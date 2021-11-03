#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <math.h>

#include "hashpipe.h"
#include "hpguppi_blade_databuf.h"
#include "hpguppi_blade.h"

#include <cuda_runtime.h>

#define ELAPSED_S(start,stop) \
  ((int64_t)stop.tv_sec-start.tv_sec)

#define ELAPSED_NS(start,stop) \
  (ELAPSED_S(start,stop)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

static int init(hashpipe_thread_args_t *args){
  cudaSetDevice(args->instance_id);
  return 0;
}

static void *run(hashpipe_thread_args_t *args)
{
  hpguppi_input_databuf_t *indb  = (hpguppi_input_databuf_t *)args->ibuf;
  hpguppi_blade_output_databuf_t *outdb = (hpguppi_blade_output_databuf_t *)args->obuf;
  char * databuf_header;

  hashpipe_status_t st = args->st;
  const char* status_key = args->thread_desc->skey;
  const char* thread_name = args->thread_desc->name;

  int curblock_in=0;
  int curblock_out=0;

  /* Sundry flags */
  int hpguppi_databuf_wait_rv, status_state=0, update_status=1;

  /* Timestamp variables */
  struct timespec ts_status_update = {0}, ts_now = {0};
  const uint64_t status_update_period_ns = 1e9;
  
  uint64_t fill_to_free_elapsed_ns;
  uint64_t fill_to_free_moving_sum_ns = 0;
  uint64_t fill_to_free_block_ns[N_INPUT_BLOCKS] = {0};
  struct timespec ts_free_input = {0};

  /* BLADE variables */
  const size_t batch_size = 2;

  struct timespec *ts_blocks_recvd = (struct timespec *)malloc(batch_size * sizeof(struct timespec));
  
  size_t inputs_batched = 0;
  size_t outputs_batched = 0;
  int *batched_input_buffer_indices = (int *)malloc(batch_size * sizeof(int));
  void **blade_input_buffers = (void**)malloc(batch_size * sizeof(void*));
  int *batched_output_buffer_indices = (int *)malloc(batch_size * sizeof(int));
  void **blade_output_buffers = (void**)malloc(batch_size * sizeof(void*));
  
  module_t mod = blade_init(batch_size);

  const size_t beamformer_input_dim_NANTS = get_input_dim_NANTS(mod);
  int32_t input_buffer_dim_NANTS, prev_flagged_NANTS;
  const size_t beamformer_input_dim_NCHANS = get_input_dim_NCHANS(mod);
  int32_t input_buffer_dim_NCHAN, prev_flagged_NCHAN;
  const size_t beamformer_input_dim_NTIME = get_input_dim_NTIME(mod);
  int32_t input_buffer_dim_NTIME, prev_flagged_NTIME;
  const size_t beamformer_input_dim_NPOLS = get_input_dim_NPOLS(mod);
  int32_t input_buffer_dim_NPOLS, prev_flagged_NPOLS;

  while (run_threads())
  {
    do {
      hpguppi_databuf_wait_rv = hpguppi_input_databuf_wait_filled(
          indb, curblock_in);
      
      clock_gettime(CLOCK_MONOTONIC, &ts_now);
      memcpy(ts_blocks_recvd+inputs_batched, &ts_now, sizeof(struct timespec));

      update_status = ELAPSED_NS(ts_status_update, ts_now) > status_update_period_ns;

      if(hpguppi_databuf_wait_rv == HASHPIPE_TIMEOUT && !update_status) {
        // No, continue receiving
        continue;
      }

      // We perform some status buffer updates every second
      if(update_status) {
        hashpipe_status_lock_safe(&st);
        {
              hputr4(st.buf, "BMBLKMS",
                round((double)fill_to_free_moving_sum_ns / N_INPUT_BLOCKS) / 1e6);
        }
        hashpipe_status_unlock_safe(&st);
      }

      // Set status field to "status_state" if we are not getting packets
      if(hpguppi_databuf_wait_rv == HASHPIPE_TIMEOUT && run_threads() && !status_state) {
        hashpipe_status_lock_safe(&st);
        {
          hputs(st.buf, status_key, "inblocked");
        }
        hashpipe_status_unlock_safe(&st);
        status_state=1;
      }
      else if (hpguppi_databuf_wait_rv != HASHPIPE_TIMEOUT && hpguppi_databuf_wait_rv != HASHPIPE_OK) {
        hashpipe_error(thread_name, "error status_state for input buffer, rv: %i", hpguppi_databuf_wait_rv);
        pthread_exit(NULL);
      }

      // Will exit if thread has been cancelled
      pthread_testcancel();
    } while (hpguppi_databuf_wait_rv != HASHPIPE_OK && run_threads()); // end wait for data loop

    // collect new input_buffer until there are `batch_size` input_buffers to process
    
    databuf_header = hpguppi_databuf_header(indb, curblock_in);
    hgeti4(databuf_header, "NANTS", &input_buffer_dim_NANTS);
    hgeti4(databuf_header, "NCHAN", &input_buffer_dim_NCHAN);
    hgeti4(databuf_header, "PIPERBLK", &input_buffer_dim_NTIME);
    hgeti4(databuf_header, "NPOL", &input_buffer_dim_NPOLS);

    if(input_buffer_dim_NANTS != beamformer_input_dim_NANTS && prev_flagged_NANTS != input_buffer_dim_NANTS){
      prev_flagged_NANTS = input_buffer_dim_NANTS;
      hashpipe_error(thread_name, "Incoming data_buffer has NANTS %lu != %lu. Ignored.", input_buffer_dim_NANTS, beamformer_input_dim_NANTS);
      hpguppi_input_databuf_set_free(indb, curblock_in);
    }
    else if(input_buffer_dim_NCHAN != beamformer_input_dim_NCHANS && prev_flagged_NCHAN != input_buffer_dim_NCHAN){
      prev_flagged_NCHAN = input_buffer_dim_NCHAN;
      hashpipe_error(thread_name, "Incoming data_buffer has NCHANS %lu != %lu. Ignored.", input_buffer_dim_NCHAN, beamformer_input_dim_NCHANS);
      hpguppi_input_databuf_set_free(indb, curblock_in);
    }
    else if(input_buffer_dim_NTIME != beamformer_input_dim_NTIME && prev_flagged_NTIME != input_buffer_dim_NTIME){
      prev_flagged_NTIME = input_buffer_dim_NTIME;
      hashpipe_error(thread_name, "Incoming data_buffer has NTIME %lu != %lu. Ignored.", input_buffer_dim_NTIME, beamformer_input_dim_NTIME);
      hpguppi_input_databuf_set_free(indb, curblock_in);
    }
    else if(input_buffer_dim_NPOLS != beamformer_input_dim_NPOLS && prev_flagged_NPOLS != input_buffer_dim_NPOLS){
      prev_flagged_NPOLS = input_buffer_dim_NPOLS;
      hashpipe_error(thread_name, "Incoming data_buffer has NPOLS %lu != %lu. Ignored.", input_buffer_dim_NPOLS, beamformer_input_dim_NPOLS);
      hpguppi_input_databuf_set_free(indb, curblock_in);
    }
    else{
      prev_flagged_NANTS = input_buffer_dim_NANTS;
      prev_flagged_NCHAN = input_buffer_dim_NCHAN;
      prev_flagged_NTIME = input_buffer_dim_NTIME;
      prev_flagged_NPOLS = input_buffer_dim_NPOLS;
      blade_input_buffers[inputs_batched] = hpguppi_databuf_data(indb, curblock_in);
      batched_input_buffer_indices[inputs_batched] = curblock_in;
      // hashpipe_info(thread_name, "batched block #%d as input #%d.", curblock_in, inputs_batched);
      inputs_batched += 1;
    }
    curblock_in  = (curblock_in + 1) % indb->header.n_block;

    if(inputs_batched < batch_size){
      continue;
    }
    
    // wait for `batch_size` output_buffers to be free

    // waiting for output
    while(outputs_batched < batch_size){
      while ((hpguppi_databuf_wait_rv=hpguppi_blade_output_databuf_wait_free(outdb, curblock_out)) !=
          HASHPIPE_OK)
      {
        if (hpguppi_databuf_wait_rv == HASHPIPE_TIMEOUT && status_state != 2)
        {
          hashpipe_status_lock_safe(&st);
          hputs(st.buf, status_key, "outblocked");
          hashpipe_status_unlock_safe(&st);
          status_state = 2;
          continue;
        }
        else
        {
          hashpipe_error(thread_name, "error waiting for output buffer, rv: %i", hpguppi_databuf_wait_rv);
          pthread_exit(NULL);
          break;
        }
      }
      // collect free output_buffer 
      if(hpguppi_databuf_wait_rv == HASHPIPE_OK){
        blade_output_buffers[outputs_batched] = hpguppi_blade_databuf_data(outdb, curblock_out);
        batched_output_buffer_indices[outputs_batched] = curblock_out;
        // hashpipe_info(thread_name, "batched block #%d as output #%d.", curblock_out, outputs_batched);
        outputs_batched += 1;
        curblock_out = (curblock_out + 1) % outdb->header.n_block;
      }
      else {
        // presume an error occurred while status_state for a free output buffer
        break;
      }
    }

    if(status_state != 3){
      hashpipe_status_lock_safe(&st);
      hputs(st.buf, status_key, "beamforming");
      hashpipe_status_unlock_safe(&st);
      status_state = 3;
    }

    blade_process(mod, blade_input_buffers, blade_output_buffers);

    for (size_t b = 0; b < batch_size; b++)
    {
      // copy across the header
      memcpy(hpguppi_blade_databuf_header(outdb, batched_output_buffer_indices[b]), 
            hpguppi_databuf_header(indb, batched_input_buffer_indices[b]), 
            HASHPIPE_STATUS_TOTAL_SIZE);	
      
      //TODO upate output_buffer headers to reflect that they contain beams
      
      hpguppi_input_databuf_set_free(indb, batched_input_buffer_indices[b]);
      hpguppi_blade_output_databuf_set_filled(outdb, batched_output_buffer_indices[b]);
      inputs_batched --;
      outputs_batched --;

      // Update moving sum (for moving average)
      clock_gettime(CLOCK_MONOTONIC, &ts_free_input);
      fill_to_free_elapsed_ns = ELAPSED_NS(ts_blocks_recvd[b], ts_free_input);
      // Add new value, subtract old value
      fill_to_free_moving_sum_ns +=
          fill_to_free_elapsed_ns - fill_to_free_block_ns[batched_input_buffer_indices[b]];
      // Store new value
      fill_to_free_block_ns[batched_input_buffer_indices[b]] = fill_to_free_elapsed_ns;
    }
  }

  hashpipe_info(thread_name, "returning");
  blade_deinit(mod);
  return NULL;
}

static hashpipe_thread_desc_t blade_beamformer_thread = {
  name: "hpguppi_ata_blade_beamformer_thread",
  skey: "BEAMSTAT",
  init: init,
  run: run,
  ibuf_desc: {hpguppi_input_databuf_create},
  obuf_desc: {hpguppi_blade_output_databuf_create}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(& blade_beamformer_thread);
}
