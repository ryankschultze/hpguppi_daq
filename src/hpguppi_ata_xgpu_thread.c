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

#include "xgpu.h"
#include <math.h>

#ifndef DEBUG_RAWSPEC_CALLBACKS
#define DEBUG_RAWSPEC_CALLBACKS (0)
#endif

#define ELAPSED_S(start,stop) \
  ((int64_t)stop.tv_sec-start.tv_sec)

#define ELAPSED_NS(start,stop) \
  (ELAPSED_S(start,stop)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

static void *run(hashpipe_thread_args_t *args)
{
  // Local aliases to shorten access to args fields
  // Our output buffer happens to be a hpguppi_input_databuf
  hpguppi_input_xgpu_databuf_t *dbin = (hpguppi_input_xgpu_databuf_t *)args->ibuf;
  hpguppi_output_xgpu_databuf_t *dbout = (hpguppi_output_xgpu_databuf_t *)args->obuf;
  hashpipe_status_t *st = &args->st;
  const char *thread_name = args->thread_desc->name;
  const char *status_key = args->thread_desc->skey;

  // Used to calculate moving average of fill-to-free times for input blocks
  uint64_t fill_to_free_elapsed_ns;
  uint64_t fill_to_free_moving_sum_ns = 0;
  uint64_t fill_to_free_block_ns[N_INPUT_BLOCKS] = {0};
  struct timespec ts_now = {0}, ts_stop_recv = {0}, ts_updated_status_info = {0};
  const uint64_t status_info_refresh_period_ns = 200*1000*1000;
  uint64_t status_info_refresh_elapsed_ns;

  //Some xGPU stuff
  XGPUInfo xgpu_info;
  int xgpu_error = 0;
  int cudaDeviceId = args->instance_id;
  hashpipe_status_lock_safe(st);
  {
    hgeti4(st->buf, "CUDADEV", &cudaDeviceId);
    hputi4(st->buf, "CUDADEV", cudaDeviceId);
  }
  hashpipe_status_unlock_safe(st);

  unsigned int ant_nchan, npol;
  uint64_t timesample_perblock;
  double corr_integration_time, tbin;
  unsigned int blocks_in_integration = 0;
  unsigned int integration_block_count = 0;
  // Get sizing info from library
  xgpuInfo(&xgpu_info);

  printf("Using xgpu: %u stations with %u channels and integration length %u\n",
         xgpu_info.nstation, xgpu_info.nfrequency, xgpu_info.ntime);

  XGPUContext xgpu_context;
  // register all blocks as one region, and use input_offset
  xgpu_context.array_h = (ComplexInput *)(dbin->block[0].data); //input; this will stop xGPU from allocating input buffer
  xgpu_context.array_len = (dbin->header.n_block*sizeof(hpguppi_input_xgpu_block_t) - BLOCK_HDR_SIZE)/sizeof(ComplexInput);//xgpu_info.vecLength;
  xgpu_context.matrix_h = (Complex *)(dbout->block[0].data); //output; direct xGPU to use ouptut databufs
  xgpu_context.matrix_len = (dbout->header.n_block*sizeof_hpguppi_output_xgpu_block_t() - BLOCK_HDR_SIZE)/sizeof(Complex);//xgpu_info.triLength;

  xgpu_error = xgpuInit(&xgpu_context, cudaDeviceId);
  if (xgpu_error)
  {
    fprintf(stderr, "xgpuInit returned error code %d\n", xgpu_error); //XXX do something
  }

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
  // uint64_t blk_start_pktidx, obs_stop_pktidx;
  int first_block = 1, obsdone = 1, obsdone_prev = 1;
  float uvh5_nsamples; // calculated for uvh5 benefit
  uint64_t blk_start_pktidx = 0, obs_start_pktidx = 0;

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
    rv = hpguppi_input_xgpu_databuf_wait_filled(dbin, curblock_in);
    clock_gettime(CLOCK_MONOTONIC, &ts_now);
    status_info_refresh_elapsed_ns = ELAPSED_NS(ts_updated_status_info, ts_now);
    if(status_info_refresh_elapsed_ns > status_info_refresh_period_ns) {
      memcpy(&ts_updated_status_info, &ts_now, sizeof(struct timespec));

      hashpipe_status_lock_safe(st);
        hputr4(st->buf, "XGPBLKMS",
                round((double)fill_to_free_moving_sum_ns / dbin->header.n_block) / 1e6);
      hashpipe_status_unlock_safe(st);
    }
    if (rv != HASHPIPE_OK)
      continue;
    memcpy(&ts_stop_recv, &ts_now, sizeof(struct timespec));

    databuf_ptr = hpguppi_xgpu_databuf_header(dbin, curblock_in);
    hgetu8(databuf_ptr, "PIPERBLK", &timesample_perblock);// coincident
    hgetu4(databuf_ptr, "NCHAN", &ant_nchan);
    hgetu4(databuf_ptr, "NPOL", &npol);

    if (xgpu_info.npol != npol ||
        //xgpu_info.nstation != nstation ||
        xgpu_info.nfrequency != ant_nchan ||
        xgpu_info.ntime != timesample_perblock)
    {
      hashpipe_error(thread_name, "Correlating %u stations with %u channels and integration length %u\n\tObservation has %u channels and %u TimeSamplesPerBlock (PIPERBLK)\n",
                      xgpu_info.nstation, xgpu_info.nfrequency, xgpu_info.ntime,
                      ant_nchan, timesample_perblock);

      /* Mark as free */
      hpguppi_input_xgpu_databuf_set_free(dbin, curblock_in);

      /* Go to next block */
      curblock_in = (curblock_in + 1) % dbin->header.n_block;

      continue;
    }

    hgetu8(databuf_ptr, "BLKSTART", &blk_start_pktidx);
    hgetu8(databuf_ptr, "PKTSTART", &obs_start_pktidx);
    if(blk_start_pktidx < obs_start_pktidx) {
      do {
        rv = hpguppi_output_xgpu_databuf_wait_free(dbout, curblock_out);
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
      memcpy(hpguppi_xgpu_output_databuf_header(dbout, curblock_out),
        databuf_ptr,
        HASHPIPE_STATUS_TOTAL_SIZE
      );
      hpguppi_output_xgpu_databuf_set_filled(dbout, curblock_out);
      curblock_out = (curblock_out + 1) % dbout->header.n_block;

      /* Mark as free */
      hpguppi_input_xgpu_databuf_set_free(dbin, curblock_in);

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
      // TODO if start obs clear the integration buffer
      first_block = 0;

      hashpipe_info(thread_name, "Observation start idx: %lu, first block start idx: %lu", obs_start_pktidx, blk_start_pktidx);

      hgetr8(databuf_ptr, "XTIMEINT", &corr_integration_time);
      hgetr8(databuf_ptr, "TBIN", &tbin);

      integration_block_count = 0;
      blocks_in_integration = (unsigned int) (corr_integration_time / (tbin * timesample_perblock) + 0.5);
      if(blocks_in_integration == 0) {
        blocks_in_integration = 1;
      }

      hashpipe_status_lock_safe(st);
        hgetu8(st->buf, "NDROP", &ndrop_integration_start);
      hashpipe_status_unlock_safe(st);

      hashpipe_status_lock_safe(st);
      hputr8(st->buf, "XTIMEINT", blocks_in_integration * timesample_perblock * tbin);
      hashpipe_status_unlock_safe(st);
      hputr8(databuf_ptr, "XTIMEINT", blocks_in_integration * timesample_perblock * tbin);

      hashpipe_info(thread_name, "Clearing integration every %u block(s).", blocks_in_integration);
    }
    hget_obsdone(st, &obsdone);
    first_block = !obsdone_prev && obsdone; // tag end of observation
    obsdone_prev = obsdone;

    do {
      rv = hpguppi_output_xgpu_databuf_wait_free(dbout, curblock_out);
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

    if(status_index !=  2)
    {
      status_index = 2;
      hashpipe_status_lock_safe(st);
      hputs(st->buf, status_key, "processing");
      hashpipe_status_unlock_safe(st);
    }

    // Correlate
    xgpu_context.input_offset = (curblock_in*sizeof(hpguppi_input_xgpu_block_t))/sizeof(ComplexInput);
    if(integration_block_count == 0) {
      // Clear the previous integration
      xgpuClearDeviceIntegrationBuffer(&xgpu_context);
      xgpu_context.output_offset = (curblock_out*sizeof_hpguppi_output_xgpu_block_t())/sizeof(Complex);
    }

    xgpu_error = xgpuCudaXengine(&xgpu_context, ++integration_block_count == blocks_in_integration ? SYNCOP_DUMP : SYNCOP_SYNC_TRANSFER); //XXX figure out flags (85 Gbps with SYNCOP_DUMP)
    if (xgpu_error)
      fprintf(stderr, "xgpuCudaXengine returned error code %d\n", xgpu_error);

    if (integration_block_count == blocks_in_integration) {
      integration_block_count = 0;
      // Correlation done. Reorder the matrix and dump to disk
      xgpuReorderMatrix(xgpu_context.matrix_h + xgpu_context.output_offset);

      hashpipe_status_lock_safe(st);
        hgetu8(st->buf, "NDROP", &ndrop_integration_end);
      hashpipe_status_unlock_safe(st);
      hputu8(databuf_ptr, "NDROP", ndrop_integration_end - ndrop_integration_start);

      uvh5_nsamples = 1.0 - ((float)(ndrop_integration_end - ndrop_integration_start))/(blocks_in_integration*timesample_perblock);
      hputr4(databuf_ptr, "NSAMPLES", uvh5_nsamples);
      hputi4(databuf_ptr, "X_TRILEN", xgpu_info.triLength); // Communicate to the downstream thread

      /* Mark as filled */
      memcpy(hpguppi_xgpu_output_databuf_header(dbout, curblock_out),
        databuf_ptr,
        HASHPIPE_STATUS_TOTAL_SIZE
      );
      hpguppi_output_xgpu_databuf_set_filled(dbout, curblock_out);
      curblock_out = (curblock_out + 1) % dbout->header.n_block;

      ndrop_integration_start = ndrop_integration_end;
    }

    /* Mark as free */
    hpguppi_input_xgpu_databuf_set_free(dbin, curblock_in);

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

static hashpipe_thread_desc_t xgpu_thread = {
  name : "hpguppi_ata_xgpu_thread",
  skey : "XGPUSTAT",
  init : NULL,
  run : run,
  ibuf_desc : {hpguppi_input_xgpu_databuf_create},
  obuf_desc : {hpguppi_output_xgpu_databuf_create}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&xgpu_thread);
}

// vi: set ts=8 sw=2 noet :