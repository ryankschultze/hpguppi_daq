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

#include "xgpu_info.h"

#include "cube/cube.h"
#include "xgpu.h"
#include <math.h>

#ifndef DEBUG_RAWSPEC_CALLBACKS
#define DEBUG_RAWSPEC_CALLBACKS (0)
#endif

#define ELAPSED_S(start,stop) \
  ((int64_t)stop.tv_sec-start.tv_sec)

#define ELAPSED_NS(start,stop) \
  (ELAPSED_S(start,stop)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

static ssize_t write_all(int fd, const void *buf, size_t bytes_to_write)
{
  size_t bytes_remaining = bytes_to_write;
  ssize_t bytes_written = 0;
  while (bytes_remaining != 0)
  {
    bytes_written = write(fd, buf, bytes_remaining);
    if (bytes_written == -1)
    {
      // Error!
      return -1;
    }
    bytes_remaining -= bytes_written;
    buf += bytes_written;
  }
  // All done!
  return bytes_to_write;
}

static int safe_close(int *pfd)
{
  if (pfd == NULL)
    return 0;
  fsync(*pfd);
  return close(*pfd);
}

int writeCorrelatorOutput(FILE *fname, Complex *buf, size_t bytes_to_write)
{
  fwrite(buf, 1, bytes_to_write, fname);
  return 1;
}

// Extend xGPU's data to include 4 redundant antennas
// Because xGPU fails with only 12
// xGPU format is:
//
//    [Slowest ---> Fastest]
//    Time        [0 ... PIPERBLK*PKTNTIME]
//    Channel     [0 ... NSTRM*PKTNCHAN]
//    FENG        [0 ... NANT]
//    POL         [0 ... NPOL]
//
//    We're injecting 4 FENGs into destination
void patch_ants(void *dest, void *src, int nants, int nants_red,
                int npols, int nchans, int ntime)
{
  char *dest_ptr = (char *)dest;
  char *src_ptr = (char *)src;

  size_t ncpy = nants * npols;
  size_t nred = nants_red * npols;

  for (size_t itime = 0; itime < ntime; itime++)
  {
    for (size_t ichan = 0; ichan < nchans; ichan++)
    {
      memcpy(dest_ptr, src_ptr, ncpy);

      dest_ptr += (ncpy + nred);
      src_ptr += ncpy;
    }
  }
}

static void *run(hashpipe_thread_args_t *args)
{
  // Local aliases to shorten access to args fields
  // Our output buffer happens to be a hpguppi_input_databuf
  hpguppi_input_xgpu_databuf_t *db = (hpguppi_input_xgpu_databuf_t *)args->ibuf;
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

  unsigned int nstation;
  uint64_t corr_time_integrations;
  unsigned int integration_clear_counter = 0;
  // Get sizing info from library
  xgpuInfo(&xgpu_info);

  printf("Using xgpu: %u stations with %u channels and integration length %u\n",
         xgpu_info.nstation, xgpu_info.nfrequency, xgpu_info.ntime);

  XGPUContext xgpu_context;
  // register all blocks as one region, and use input_offset
  xgpu_context.array_h = (ComplexInput *)(db->block[0].hdr); //input; this will stop xGPU from allocating input buffer
  xgpu_context.array_len = db->header.n_block*sizeof(hpguppi_input_xgpu_block_t)/sizeof(ComplexInput);//xgpu_info.vecLength;
  xgpu_context.matrix_h = NULL;           //output; xGPU will allocate memory and take care of this internally

  xgpu_error = xgpuInit(&xgpu_context, cudaDeviceId);
  if (xgpu_error)
  {
    fprintf(stderr, "xgpuInit returned error code %d\n", xgpu_error); //XXX do something
  }

  Complex *cuda_matrix_h = xgpu_context.matrix_h;
  ///////////////////////////

  /* Read in general parameters */
  struct hpguppi_params gp;
  struct psrfits pf;
  pf.sub.dat_freqs = NULL;
  pf.sub.dat_weights = NULL;
  pf.sub.dat_offsets = NULL;
  pf.sub.dat_scales = NULL;
  pthread_cleanup_push((void *)hpguppi_free_psrfits, &pf);

  /* Init output file descriptor (-1 means no file open) */
  static int fdraw = -1;
  pthread_cleanup_push((void *)safe_close, &fdraw);

  /* Set I/O priority class for this thread to "real time" */
  if (ioprio_set(IOPRIO_WHO_PROCESS, 0, IOPRIO_PRIO_VALUE(IOPRIO_CLASS_RT, 7)))
  {
    hashpipe_error(thread_name, "ioprio_set IOPRIO_CLASS_RT");
  }

  /* Loop */
  uint64_t pktidx = 0, pktstart = 0, pktstop = 0;
  int len = 0;
  int curblock = 0;
  int block_count = 0, blocks_per_file = 128, filenum = 0;
  int got_packet_0 = 0, first = 1;
  char *ptr; //, *hend;
  int open_flags = 0;
  int directio = 0;
  char fname[256];
  unsigned char base_filename_stem_start;
  int rv = 0;
  int status_index = -1;
  //char* tmp_data = malloc(xgpu_info.vecLength*sizeof(char)); //XXX

  int mjd_d, mjd_s;
  double mjd_fs;

  /* Heartbeat variables */
  time_t curtime = 0;
  char timestr[32] = {0};
  enum run_states state = IDLE;

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
    rv = hpguppi_input_xgpu_databuf_wait_filled(db, curblock);
    clock_gettime(CLOCK_MONOTONIC, &ts_now);
    status_info_refresh_elapsed_ns = ELAPSED_NS(ts_updated_status_info, ts_now);
    if(status_info_refresh_elapsed_ns > status_info_refresh_period_ns) {
      memcpy(&ts_updated_status_info, &ts_now, sizeof(struct timespec));

      time(&curtime);//time stores seconds since epoch
      ctime_r(&curtime, timestr);
      timestr[strlen(timestr)-1] = '\0'; // Chop off trailing newline

      hashpipe_status_lock_safe(st);
        hputr4(st->buf, "XGPBLKMS",
                round((double)fill_to_free_moving_sum_ns / db->header.n_block) / 1e6);
        hputs(st->buf, "DAQSTATE", state == IDLE  ? "idling" :\
                             state == ARMED ? "armed"  :\
                             "recording");
        hputs(st->buf, "DAQPULSE", timestr);
      hashpipe_status_unlock_safe(st);
    }
    if (rv != HASHPIPE_OK)
      continue;
    memcpy(&ts_stop_recv, &ts_now, sizeof(struct timespec));

    /* Read param struct for this block */
    ptr = hpguppi_xgpu_databuf_header(db, curblock);
    if (first)
    {
      hpguppi_read_obs_params(ptr, &gp, &pf);
      hashpipe_info(thread_name, "First block");
      first = 0;
    }
    else
    {
      hpguppi_read_subint_params(ptr, &gp, &pf);

      // mjd_d = pf.hdr.start_day;
      // mjd_s = 0;
      // mjd_fs = pf.hdr.start_sec;

      hgeti4(ptr, "STT_IMJD", &mjd_d);
      hgeti4(ptr, "STT_SMJD", &mjd_s);
      hgetr8(ptr, "STT_OFFS", &mjd_fs);

      if (pf.hdr.start_day != mjd_d || pf.hdr.start_sec != mjd_s + mjd_fs)
      { // Observation timestamp has changed. Start new stem.
        if (fdraw != -1)
        {
          fprintf(stderr, "STT_MJD value changed. Starting a new file stem. (Previous stem had %d files)\n", filenum);
          close(fdraw);
          // Reset fdraw, got_packet_0, filenum, block_count
          fdraw = -1;
          got_packet_0 = 0;
          filenum = 0;
          block_count = 0;
        }
      }
    }

    /* Read pktidx, pktstart, pktstop from header */
    hgetu8(ptr, "PKTIDX", &pktidx);
    hgetu8(ptr, "PKTSTART", &pktstart);
    hgetu8(ptr, "PKTSTOP", &pktstop);

    // If packet idx is NOT within start/stop range
    if (pktidx < pktstart || pktstop <= pktidx)
    {
      state = pktidx < pktstart ? ARMED : IDLE;
      // If file open, close it
      if (fdraw != -1)
      {
        // Close file
        close(fdraw);
        // Reset fdraw, got_packet_0, filenum, block_count
        fdraw = -1;
        got_packet_0 = 0;
        filenum = 0;
        block_count = 0;
        hput_obsdone(st, 1);

        // Print end of recording conditions
        hashpipe_info(thread_name, "recording stopped: "
                                   "pktstart %lu pktstop %lu pktidx %lu",
                      pktstart, pktstop, pktidx);
        fprintf(stderr, "xgpu thread: status buffer follows\n");

        hgetu8(ptr, "PKTIDX", &pktidx);
        fprintf(stderr, "hgetu8(PKTIDX) buffered: %lu\n", pktidx);
        fprintf(stderr, "buffered:\n%s", ptr);

        fprintf(stderr, "\nlive ptr below:\n");

        hgetu8(ptr, "PKTIDX", &pktidx);
        fprintf(stderr, "hgetu8(PKTIDX): %lu\n", pktidx);
        fprintf(stderr, "%s", ptr);
      }
      else
      {
        hashpipe_info(thread_name, "Block not recorded: "
                                   "pktstart %lu pktstop %lu pktidx %lu",
                      pktstart, pktstop, pktidx);
      }
      /* Mark as free */
      hpguppi_input_xgpu_databuf_set_free(db, curblock);

      /* Go to next block */
      curblock = (curblock + 1) % db->header.n_block;

      continue;
    }
    state = RECORD;

    // Wait for packet 0 before starting write
    // "packet 0" is the first packet/block of the new recording,
    // it is not necessarily pktidx == 0.
    if (got_packet_0 == 0)
    { // && gp.stt_valid==1) {
      got_packet_0 = 1;
      hput_obsdone(st, 1);

      hpguppi_read_obs_params(ptr, &gp, &pf);
      corr_time_integrations = gp.packets_per_block;
      hgetu4(ptr, "NANTS", &nstation);
      hgetu8(ptr, "XTIMEINT", &corr_time_integrations);
      hputu8(ptr, "XTIMEINT", corr_time_integrations);
      integration_clear_counter = corr_time_integrations/gp.packets_per_block; // trigger a clearing of the integration
      hashpipe_info(thread_name, "Clearing integration every %u block(s).", integration_clear_counter);

      if (xgpu_info.npol != pf.hdr.npol ||
          //xgpu_info.nstation != nstation ||
          xgpu_info.nfrequency != pf.hdr.nchan / nstation ||
          xgpu_info.ntime != gp.packets_per_block)
      {
        hashpipe_error(thread_name, "Correlating %u stations with %u channels and integration length %u\n\tObservation has %u stations with %u channels and TimeSamplesPerBlock (PIPERBLK) %u\n",
                       xgpu_info.nstation, xgpu_info.nfrequency, xgpu_info.ntime,
                       nstation, pf.hdr.nchan / nstation, gp.packets_per_block);
        got_packet_0 = 0;

        /* Mark as free */
        hpguppi_input_xgpu_databuf_set_free(db, curblock);

        /* Go to next block */
        curblock = (curblock + 1) % db->header.n_block;

        continue;
      }

      directio = hpguppi_read_directio_mode(ptr);
      sprintf(fname, "%s.%04d.xgpu", pf.basefilename, filenum);
      fprintf(stderr, "Opening first file '%s' (directio=%d)\n", fname, directio);

      // finds last '/'
      base_filename_stem_start = strlen(pf.basefilename);
      while(base_filename_stem_start > 0 && pf.basefilename[base_filename_stem_start-1] != '/'){
        base_filename_stem_start--;
      }

      hashpipe_status_lock_safe(st);
        hputs(st->buf, "OBSSTEM", pf.basefilename+base_filename_stem_start);
      hashpipe_status_unlock_safe(st);

      // Create the output directory if needed
      char datadir[1024];
      strncpy(datadir, pf.basefilename, 1023);
      char *last_slash = strrchr(datadir, '/');
      if (last_slash != NULL && last_slash != datadir)
      {
        *last_slash = '\0';
        printf("Using directory '%s' for output.\n", datadir);
        if (mkdir_p(datadir, 0755) == -1)
        {
          hashpipe_error(thread_name, "mkdir_p(%s)", datadir);
          break;
        }
      }
      // TODO: check for file exist.
      open_flags = O_CREAT | O_RDWR;
      if (directio)
      {
        open_flags |= O_DIRECT;
      }
      fdraw = open(fname, open_flags, 0644);
      if (fdraw == -1)
      {
        hashpipe_error(thread_name, "Error opening file.");
        pthread_exit(NULL);
      }
    }

    /* See if we need to open next file */
    if (block_count >= blocks_per_file)
    {
      close(fdraw);
      filenum++;
      char fname[256];
      sprintf(fname, "%s.%4.4d.xgpu", pf.basefilename, filenum);
      directio = hpguppi_read_directio_mode(ptr);
      open_flags = O_CREAT | O_RDWR;
      if (directio)
      {
        open_flags |= O_DIRECT;
      }
      fprintf(stderr, "Opening next file '%s' (directio=%d)\n", fname, directio);
      fdraw = open(fname, open_flags, 0644);
      if (fdraw == -1)
      {
        hashpipe_error(thread_name, "Error opening file.");
        pthread_exit(NULL);
      }
      block_count = 0;
    }

    /* If we got packet 0, write data to disk */
    if (got_packet_0)
    {
      /* Note writing status */
      if(status_index !=  1)
      {
        status_index = 1;
        hashpipe_status_lock_safe(st);
        hputs(st->buf, status_key, "processing");
        hashpipe_status_unlock_safe(st);
      }

      //memcpy(tmp_data, pf.sub.data, BLOCK_DATA_SIZE);

      // fprintf(stderr, "xGPU thread [%d]: %i %i %i %i %i\n", curblock, pf.sub.data[0],
      // 	  pf.sub.data[100], pf.sub.data[200], pf.sub.data[300],
      // 	  pf.sub.data[400]);

      xgpu_context.input_offset = (db->block[curblock].data - db->block[0].hdr)/sizeof(ComplexInput);
      //xgpu_context.array_h = (ComplexInput*) tmp_data; //XXX

      // Clear the previous integration
      if(integration_clear_counter == corr_time_integrations/gp.packets_per_block)
      {
        xgpuClearDeviceIntegrationBuffer(&xgpu_context);
        integration_clear_counter = 0;
      }

      // Correlate
      xgpu_error = xgpuCudaXengine(&xgpu_context, SYNCOP_DUMP); //XXX figure out flags (85 Gbps with SYNCOP_DUMP)
      if (xgpu_error)
        fprintf(stderr, "xgpuCudaXengine returned error code %d\n", xgpu_error);

      // Correlation done. Reorder the matrix and dump to disk
      xgpuReorderMatrix(cuda_matrix_h);
      integration_clear_counter++;

      /* Note writing status */
      if(status_index !=  2)
      {
        status_index = 2;
        hashpipe_status_lock_safe(st);
        hputs(st->buf, status_key, "writing");
        hashpipe_status_unlock_safe(st);
      }

      //len = xgpu_context.matrix_len*sizeof(Complex); //
      len = xgpu_info.triLength * sizeof(Complex);
      //int triLength_valid = pf.hdr.npol * pf.hdr.npol * nstation *
      //	  (nstation + 1) / 2 * (pf.hdr.nchan/nstation);
      //len = triLength_valid *sizeof(Complex);

      // Adjust length for any padding required for DirectIO
      if (directio)
      {
        // Round up to next multiple of 512
        len = (len + 511) & ~511;
      }

      rv = write_all(fdraw, cuda_matrix_h, len);
      if (rv != len)
      {
        char msg[100];
        perror(thread_name);
        sprintf(msg, "Error writing data (ptr=%p, len=%d, rv=%d)", ptr, len, rv);
        hashpipe_error(thread_name, msg);
      }

      if (!directio)
      {
        /* flush output */
        fsync(fdraw);
      }

      /* Increment counter */
      block_count++;
    }

    /* Mark as free */
    hpguppi_input_xgpu_databuf_set_free(db, curblock);

    // Update moving sum (for moving average)
    clock_gettime(CLOCK_MONOTONIC, &ts_now);
    fill_to_free_elapsed_ns = ELAPSED_NS(ts_stop_recv, ts_now);
    // Add new value, subtract old value
    fill_to_free_moving_sum_ns +=
        fill_to_free_elapsed_ns - fill_to_free_block_ns[curblock];
    // Store new value
    fill_to_free_block_ns[curblock] = fill_to_free_elapsed_ns;

    /* Go to next block */
    curblock = (curblock + 1) % db->header.n_block;

    /* Check for cancel */
    pthread_testcancel();
  }

  hashpipe_info(thread_name, "exiting!");
  pthread_exit(NULL);

  pthread_cleanup_pop(0); /* Closes safe_close */
  pthread_cleanup_pop(0); /* Closes hpguppi_free_psrfits */
}

static hashpipe_thread_desc_t xgpu_thread = {
  name : "hpguppi_atasnap_xgpu_thread",
  skey : "XGPUSTAT",
  init : NULL,
  run : run,
  ibuf_desc : {hpguppi_input_xgpu_databuf_create},
  obuf_desc : {NULL}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&xgpu_thread);
}

// vi: set ts=8 sw=2 noet :
