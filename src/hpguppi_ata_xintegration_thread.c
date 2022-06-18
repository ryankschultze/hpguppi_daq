/* hpguppi_rawdisk_only_thread.c
 *
 * Write databuf blocks out to disk.
 */

#define _GNU_SOURCE 1
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#include "ioprio.h"

#include "hashpipe.h"
#include "hpguppi_atasnap.h"

#include "hpguppi_xgpu_databuf.h"
#include "hpguppi_params.h"
//#include "hpguppi_pksuwl.h"
#include "hpguppi_util.h"

#include "uvh5.h"

static void *run(hashpipe_thread_args_t *args)
{
  // Local aliases to shorten access to args fields
  // Our output buffer happens to be a hpguppi_input_databuf
  hpguppi_output_xgpu_databuf_t *dbin = (hpguppi_output_xgpu_databuf_t *)args->ibuf;
  hpguppi_output_xgpu_databuf_t *dbout = (hpguppi_output_xgpu_databuf_t *)args->obuf;
  hashpipe_status_t *st = &args->st;
  const char *thread_name = args->thread_desc->name;
  const char *status_key = args->thread_desc->skey;

  // Used to calculate moving average of fill-to-free times for input blocks
  uint64_t fill_to_free_elapsed_ns;
  uint64_t fill_to_free_moving_sum_ns = 0;
  uint64_t fill_to_free_block_ns[N_XGPU_OUTPUT_BLOCKS] = {0};
  struct timespec ts_now = {0}, ts_stop_recv = {0}, ts_updated_status_info = {0};
  const uint64_t status_info_refresh_period_ns = 200*1000*1000;
  uint64_t status_info_refresh_elapsed_ns;

  unsigned int ant_nchan, npol;
  uint64_t timesample_perblock;
  double corr_integration_time, tbin;
  unsigned int blocks_in_integration = 0;
  unsigned int integration_block_count = 0;
  unsigned int gpu_integration_block_count = 0;

  /* Set I/O priority class for this thread to "real time" */
  if (ioprio_set(IOPRIO_WHO_PROCESS, 0, IOPRIO_PRIO_VALUE(IOPRIO_CLASS_RT, 7)))
  {
    hashpipe_error(thread_name, "ioprio_set IOPRIO_CLASS_RT");
  }

  /* Loop */
  int curblock_in = 0;
  int curblock_out = 0;
  char* databuf_ptr;
  int rv = 0;
  int status_index = -1;
  uint64_t ndrop_integration_start, ndrop_integration_end;
  int first_block = 1;
  float uvh5_nsamples; // calculated for uvh5 benefit
  uint64_t blk_start_pktidx = 0, obs_start_pktidx = 0, obs_stop_pktidx = 0;
  uint64_t integrated_blk_start_pktidx = 0;
  int32_t triLength;
  UVH5_CI32_t* databuf_input_complex_ptr;
  UVH5_CF64_t* databuf_output_complex_ptr;

  while (run_threads())
  {

    /* Note waiting status */
    if(status_index !=  0)
    {
      status_index = 0;
      hashpipe_status_lock_safe(st);
      hputs(st->buf, status_key, "waiting");
      hashpipe_status_unlock_safe(st);
    }

    /* Wait for buf to have data */
    rv = hpguppi_databuf_wait_filled(dbin, curblock_in);
    clock_gettime(CLOCK_MONOTONIC, &ts_now);
    status_info_refresh_elapsed_ns = ELAPSED_NS(ts_updated_status_info, ts_now);
    if(status_info_refresh_elapsed_ns > status_info_refresh_period_ns) {
      memcpy(&ts_updated_status_info, &ts_now, sizeof(struct timespec));

      hashpipe_status_lock_safe(st);
        hputr4(st->buf, "XINBLKMS",
                round((double)fill_to_free_moving_sum_ns / dbin->header.n_block) / 1e6);
      hashpipe_status_unlock_safe(st);
    }
    if (rv != HASHPIPE_OK)
      continue;
    memcpy(&ts_stop_recv, &ts_now, sizeof(struct timespec));

    databuf_ptr = hpguppi_databuf_header(dbin, curblock_in);
    hgetu8(databuf_ptr, "PIPERBLK", &timesample_perblock);// coincident
    hgetu4(databuf_ptr, "NCHAN", &ant_nchan);
    hgetu4(databuf_ptr, "NPOL", &npol);

    hgetu8(databuf_ptr, "BLKSTART", &blk_start_pktidx);
    hgetu8(databuf_ptr, "PKTSTART", &obs_start_pktidx);
    hgetu8(databuf_ptr, "PKTSTOP", &obs_stop_pktidx);
    if(blk_start_pktidx < obs_start_pktidx
      || blk_start_pktidx > obs_stop_pktidx
    ) {
      first_block = 1;
      do {
        rv = hpguppi_databuf_wait_free(dbout, curblock_out);
        if (rv == HASHPIPE_TIMEOUT) {
          if(status_index !=  1)
          {
            status_index = 1;
            hashpipe_status_lock_safe(st);
            hputs(st->buf, status_key, "blocked");
            hashpipe_status_unlock_safe(st);
          }
        }
      } while (rv != HASHPIPE_OK);

      /* Mark useless block as filled */
      memcpy(hpguppi_databuf_header(dbout, curblock_out),
        databuf_ptr,
        HASHPIPE_STATUS_TOTAL_SIZE
      );
      hpguppi_databuf_set_filled(dbout, curblock_out);
      curblock_out = (curblock_out + 1) % dbout->header.n_block;

      /* Mark as free */
      hpguppi_databuf_set_free(dbin, curblock_in);

      // Update moving sum (for moving average)
      clock_gettime(CLOCK_MONOTONIC, &ts_now);
      fill_to_free_elapsed_ns = ELAPSED_NS(ts_stop_recv, ts_now);
      // Add new value, subtract old value
      fill_to_free_moving_sum_ns +=
          fill_to_free_elapsed_ns - fill_to_free_block_ns[curblock_in];
      // Store new value
      fill_to_free_block_ns[curblock_in] = fill_to_free_elapsed_ns;

      /* Go to next block */
      curblock_in = (curblock_in + 1) % dbin->header.n_block;
      continue;
    }

    if (first_block) {
      first_block = 0;
      integration_block_count = 0;

      hgeti4(databuf_ptr, "X_TRILEN", &triLength);
      hgetr8(databuf_ptr, "XTIMEINT", &corr_integration_time);
      hgetr8(databuf_ptr, "TBIN", &tbin);
      gpu_integration_block_count = 0;
      hgetu4(databuf_ptr, "XGPUINT", &gpu_integration_block_count);

      integration_block_count = 0;
      blocks_in_integration = (unsigned int) (corr_integration_time / (tbin * timesample_perblock) + 0.5);
      if(gpu_integration_block_count > 0){
        blocks_in_integration /= gpu_integration_block_count;
      }
      else {
        gpu_integration_block_count = blocks_in_integration;
        blocks_in_integration = 1;
      }
      if(blocks_in_integration == 0) {
        blocks_in_integration = 1;
      }
      hashpipe_info(thread_name, "Clearing CPU integration every %u block(s).", blocks_in_integration);
    }

    if(status_index !=  2)
    {
      status_index = 2;
      hashpipe_status_lock_safe(st);
      hputs(st->buf, status_key, "processing");
      hashpipe_status_unlock_safe(st);
    }

    // Integrate block into output
    databuf_input_complex_ptr = (UVH5_CI32_t*) hpguppi_databuf_data(dbin, curblock_in);
    if(integration_block_count == 0) {
      /* Wait for output buffer to be free */
      do {
        rv = hpguppi_databuf_wait_free(dbout, curblock_out);
        if (rv == HASHPIPE_TIMEOUT) {
          if(status_index !=  1)
          {
            status_index = 1;
            hashpipe_status_lock_safe(st);
            hputs(st->buf, status_key, "blocked");
            hashpipe_status_unlock_safe(st);
          }
        }
      } while (rv != HASHPIPE_OK);

      databuf_output_complex_ptr = (UVH5_CF64_t*) hpguppi_databuf_data(dbout, curblock_out);
      memset(databuf_output_complex_ptr,
        0,
        hpguppi_output_xgpu_block_data_byte_size()
      );
      
      integrated_blk_start_pktidx = blk_start_pktidx;
    }

    // Integrate INT32_complex input to DOUBLE_complex output.
    for(int i = 0; i < triLength; i++){
      databuf_output_complex_ptr[i].r += (double)databuf_input_complex_ptr[i].r;
      databuf_output_complex_ptr[i].i += (double)databuf_input_complex_ptr[i].i;
    }
    integration_block_count++;

    if (integration_block_count == blocks_in_integration) {
      integration_block_count = 0;

      hashpipe_status_lock_safe(st);
        hgetu8(st->buf, "NDROP", &ndrop_integration_end);
      hashpipe_status_unlock_safe(st);

      hputu8(databuf_ptr, "NDROP", ndrop_integration_end - ndrop_integration_start);

      uvh5_nsamples = 1.0 - ((float)(ndrop_integration_end - ndrop_integration_start))/(blocks_in_integration*gpu_integration_block_count*timesample_perblock);
      hputr4(databuf_ptr, "NSAMPLES", uvh5_nsamples);
      
      /* Mark as filled */
      memcpy(hpguppi_databuf_header(dbout, curblock_out),
        databuf_ptr,
        HASHPIPE_STATUS_TOTAL_SIZE
      );
      hputu8(hpguppi_databuf_header(dbout, curblock_out), "BLKSTART", integrated_blk_start_pktidx);
      hpguppi_databuf_set_filled(dbout, curblock_out);
      curblock_out = (curblock_out + 1) % dbout->header.n_block;

      ndrop_integration_start = ndrop_integration_end;
    }

    /* Mark as free */
    hpguppi_databuf_set_free(dbin, curblock_in);

    // Update moving sum (for moving average)
    clock_gettime(CLOCK_MONOTONIC, &ts_now);
    fill_to_free_elapsed_ns = ELAPSED_NS(ts_stop_recv, ts_now);
    // Add new value, subtract old value
    fill_to_free_moving_sum_ns +=
        fill_to_free_elapsed_ns - fill_to_free_block_ns[curblock_in];
    // Store new value
    fill_to_free_block_ns[curblock_in] = fill_to_free_elapsed_ns;

    /* Go to next block */
    curblock_in = (curblock_in + 1) % dbin->header.n_block;

    /* Check for cancel */
    pthread_testcancel();
  }

  hashpipe_info(thread_name, "exiting!");
  pthread_exit(NULL);
}

static hashpipe_thread_desc_t xintegration_thread = {
  name : "hpguppi_ata_xintegration_thread",
  skey : "XINTSTAT",
  init : NULL,
  run : run,
  ibuf_desc : {hpguppi_output_xgpu_databuf_create},
  obuf_desc : {hpguppi_output_xgpu_databuf_create}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&xintegration_thread);
}

// vi: set ts=8 sw=2 noet :
