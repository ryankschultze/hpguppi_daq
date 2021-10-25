/* hpguppi_net_thread.c
 *
 * Routine to read packets from network and put them
 * into shared memory blocks.
 */

#define _GNU_SOURCE 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "hashpipe.h"

#include "hpguppi_databuf.h"
#include "hpguppi_params.h"
#include "hpguppi_udp.h"
#include "hpguppi_time.h"
#include "hpguppi_atasnap.h"

#include "ioprio.h"

// 80 character string for the BACKEND header record.
static const char BACKEND_RECORD[] =
// 0000000000111111111122222222223333333333
// 0123456789012345678901234567890123456789
  "BACKEND = 'GUPPI   '                    " \
  "                                        ";

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
  while(bytes_remaining != 0) {
    bytes_written = write(fd, buf, bytes_remaining);
    if(bytes_written == -1) {
      // Error!
      return -1;
    }
    bytes_remaining -= bytes_written;
    buf += bytes_written;
  }
  // All done!
  return bytes_to_write;
}

static int safe_close(int *pfd) {
    if (pfd==NULL) return 0;
    fsync(*pfd);
    return close(*pfd);
}

#define HPUT_DAQ_STATE(st, state)\
  hputs(st->buf, "DAQSTATE", state == IDLE  ? "idling" :\
                             state == ARMED ? "armed"  :\
                             "recording")

static void *run(hashpipe_thread_args_t * args)
{
    // Local aliases to shorten access to args fields
    // Our output buffer happens to be a hpguppi_input_databuf
    
  hpguppi_input_databuf_t *indb  = (hpguppi_input_databuf_t *)args->ibuf;

  hashpipe_status_t* st = &(args->st);
  const char* status_key = args->thread_desc->skey;
  const char* thread_name = args->thread_desc->name;

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
  if(ioprio_set(IOPRIO_WHO_PROCESS, 0, IOPRIO_PRIO_VALUE(IOPRIO_CLASS_RT, 7))) {
    hashpipe_error(thread_name, "ioprio_set IOPRIO_CLASS_RT");
  }

  int rv = 0;
  int curblock_in=0;

  char *datablock_header;
  int blocksize=0, len=0;
  int block_count=0, blocks_per_file= (int) (((uint64_t)16<<30)/BLOCK_DATA_SIZE) , filenum=0;
  int got_packet_0=0;
  char *hend;
  int open_flags = 0;
  int directio = 0;

  blocks_per_file = blocks_per_file == 0 ? 1 :  blocks_per_file;
  hashpipe_info(thread_name, "Blocks per file: %d", blocks_per_file);

  unsigned char base_filename_stem_start;

  struct mjd_t *mjd = malloc(sizeof(struct mjd_t));

  /* Misc counters, etc */
  uint64_t obs_npacket_total=0, obs_ndrop_total=0;
  uint32_t block_npacket=0, block_ndrop=0;

  uint64_t obs_start_pktidx = 0, obs_stop_pktidx = 0;
  uint64_t block_start_pktidx = 0, block_stop_pktidx = 0;
  
  char waiting=-1, flag_state_update=0;
  enum run_states state = IDLE;

  // Used to calculate moving average of fill-to-free times for input blocks
  uint64_t fill_to_free_elapsed_ns;
  uint64_t fill_to_free_moving_sum_ns = 0;
  uint64_t fill_to_free_block_ns[N_INPUT_BLOCKS] = {0};
  struct timespec ts_free_input = {0}, ts_block_recvd = {0};
  
  /* Heartbeat variables */
  time_t lasttime = 0;
  time_t curtime = 0;
  char timestr[32] = {0};
    
  uint32_t blocks_per_second = 0;

  // Reset STTVALID so that update_stt_status_keys works
  hashpipe_status_lock_safe(st);
    hputu4(st->buf, "STTVALID", 0);
  hashpipe_status_unlock_safe(st);

  /* Main loop */
  while (run_threads()) {

    /* Wait for data */
    do {
      // Heartbeat update?
      time(&curtime);//time stores seconds since epoch
      if(flag_state_update || curtime > lasttime) {// once per second
          flag_state_update = 0;
          lasttime = curtime;

          ctime_r(&curtime, timestr);
          timestr[strlen(timestr)-1] = '\0'; // Chop off trailing newline
          hashpipe_status_lock_safe(st);
          {
              hputu8(st->buf, "OBSNPKTS", obs_npacket_total);
              hputu8(st->buf, "OBSNDROP", obs_ndrop_total);
              hputu4(st->buf, "OBSBLKPS", blocks_per_second);
              hputr4(st->buf, "OBSBLKMS",
                round((double)fill_to_free_moving_sum_ns / N_INPUT_BLOCKS) / 1e6);
              hputs(st->buf, "DAQPULSE", timestr);
              HPUT_DAQ_STATE(st, state);
          }
          hashpipe_status_unlock_safe(st);
          blocks_per_second = 0;
      }

      // Waiting for input
      rv=hpguppi_input_databuf_wait_filled(indb, curblock_in);
      clock_gettime(CLOCK_MONOTONIC, &ts_block_recvd);
      if (rv == HASHPIPE_TIMEOUT)
      {
        if(waiting != 1){
          hashpipe_status_lock_safe(st);
            hputs(st->buf, status_key, "waiting");
          hashpipe_status_unlock_safe(st);
          waiting=1;
        }
      }
      else if(rv != HASHPIPE_OK)
      {
        hashpipe_error(thread_name, "error waiting for input buffer, rv: %i", rv);
        pthread_exit(NULL);
      }

    } while (rv != HASHPIPE_OK && run_threads());

    if(!run_threads()) {
      break;
    }
    
    /* Update status if needed */
    if (waiting!=0) {
        hashpipe_status_lock_safe(st);
        hputs(st->buf, status_key, "processing");
        hashpipe_status_unlock_safe(st);
        waiting=0;
    }

    datablock_header = hpguppi_databuf_header(indb, curblock_in);
    hgetu8(datablock_header, "PKTSTART", &obs_start_pktidx);
    hgetu8(datablock_header, "PKTSTOP", &obs_stop_pktidx);
    hgetu8(datablock_header, "BLKSTART", &block_start_pktidx);
    hgetu8(datablock_header, "BLKSTOP", &block_stop_pktidx);

    switch(state_from_block_start_stop(obs_start_pktidx, obs_stop_pktidx, block_start_pktidx, block_stop_pktidx)){
      case IDLE:// If should IDLE, 
        if(state != IDLE){
          if(state == RECORD){//and recording, finalise block
            // If file open, close it
            if(fdraw != -1) {
              // Close file
              close(fdraw);
              // Reset fdraw, got_packet_0, filenum, block_count
              fdraw = -1;
              got_packet_0 = 0;
              filenum = 0;
              block_count=0;

              // Print end of recording conditions
              hashpipe_info(thread_name, "recording stopped: "
                "obs_start %lu obs_stop %lu blk_start_pktidx %lu blk_stop_pktidx %lu",
                obs_start_pktidx, obs_stop_pktidx, block_start_pktidx, block_stop_pktidx);
            }
          }
          flag_state_update = 1;
          state = IDLE;
          // Update STT with state = IDLE, setting STTVALID = 0
          update_stt_status_keys(st, state, obs_start_pktidx, mjd);
        }
        break;
      case RECORD:// If should RECORD
        if (state != RECORD){
          obs_npacket_total = 0;
          obs_ndrop_total = 0;
          if(state != ARMED){// didn't arm correctly
            update_stt_status_keys(st, state, obs_start_pktidx, mjd);
            hputu4(datablock_header, "STTVALID", 1);
            hputu4(datablock_header, "STT_IMJD", mjd->stt_imjd);
            hputu4(datablock_header, "STT_SMJD", mjd->stt_smjd);
            hputr8(datablock_header, "STT_OFFS", mjd->stt_offs);
          }
          flag_state_update = 1;
          state = RECORD;
          hgeti4(datablock_header, "STTVALID", &gp.stt_valid);
          /* Get full data block size */
          hgeti4(datablock_header, "BLOCSIZE", &blocksize);
        }

        hgetu4(datablock_header, "NPKT", &block_npacket);
        hgetu4(datablock_header, "NDROP", &block_ndrop);
        obs_npacket_total += block_npacket;
        obs_ndrop_total += block_ndrop;
        break;
      case ARMED:// If should ARM,
        if(state != ARMED){
          flag_state_update = 1;
          state = ARMED;
          update_stt_status_keys(st, state, obs_start_pktidx, mjd);
        }
      default:
        break;
    }
    
    if(state == RECORD){
      /*** RAW Disk write out BEGIN */
      hpguppi_read_subint_params(datablock_header, &gp, &pf);
      
      if (fdraw != -1 && // currently writing to file
          (pf.hdr.start_day != mjd->stt_imjd 
          || pf.hdr.start_sec != mjd->stt_smjd + mjd->stt_offs)// and Observation timestamp has changed
          ){// Start new stem.
          fprintf(stderr, "STT_MJD value changed. Starting a new file stem.\n");
          close(fdraw);
          // Reset fdraw, got_packet_0, filenum, block_count
          fdraw = -1;
          got_packet_0 = 0;
          filenum = 0;
          block_count=0;

          // Copy over current STT values
          mjd->stt_imjd = pf.hdr.start_day;
          mjd->stt_smjd = (int) pf.hdr.start_sec; // floor
          mjd->stt_offs = pf.hdr.start_sec - mjd->stt_smjd;
      }

      /* Set up data ptr for quant routines */
      pf.sub.data = (unsigned char *)hpguppi_databuf_data(indb, curblock_in);

      // Wait for packet 0 before starting write
      // "packet 0" is the first packet/block of the new recording,
      // it is not necessarily pktidx == 0.
      if (got_packet_0==0 && gp.stt_valid==1) {
        got_packet_0 = 1;
        hpguppi_read_obs_params(datablock_header, &gp, &pf);
        directio = hpguppi_read_directio_mode(datablock_header);
        char fname[256];
        sprintf(fname, "%s.%04d.raw", pf.basefilename, filenum);
        fprintf(stderr, "Opening first raw file '%s' (directio=%d)\n", fname, directio);
        
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
        if (last_slash!=NULL && last_slash!=datadir) {
          *last_slash = '\0';
          printf("Using directory '%s' for output.\n", datadir);
          if(mkdir_p(datadir, 0755) == -1) {
            hashpipe_error(thread_name, "mkdir_p(%s)", datadir);
            break;
          }
        }
        // TODO: check for file exist.
        open_flags = O_CREAT|O_WRONLY;//|O_SYNC;
        if(directio) {
          open_flags |= O_DIRECT;
        }
        fdraw = open(fname, open_flags, 0644);
        if (fdraw==-1) {
          hashpipe_error(thread_name, "Error opening file.");
          pthread_exit(NULL);
        }

      }

      /* See if we need to open next file */
      if (block_count >= blocks_per_file) {
        close(fdraw);
        filenum++;
        char fname[256];
        sprintf(fname, "%s.%4.4d.raw", pf.basefilename, filenum);
        directio = hpguppi_read_directio_mode(datablock_header);
        open_flags = O_CREAT|O_WRONLY;//|O_SYNC;
        if(directio) {
          open_flags |= O_DIRECT;
        }
        fprintf(stderr, "Opening next raw file '%s' (directio=%d)\n", fname, directio);
        fdraw = open(fname, open_flags, 0644);
        if (fdraw==-1) {
          hashpipe_error(thread_name, "Error opening file.");
          pthread_exit(NULL);
        }
        block_count=0;
      }

      /* If we got packet 0, write data to disk */
      if (got_packet_0) {
        if(waiting != -1){
          /* Note writing status */
          waiting = -1;
          hashpipe_status_lock_safe(st);
          hputs(st->buf, status_key, "writing");
          hashpipe_status_unlock_safe(st);
        }

        /* Write header to file */
        hend = ksearch(datablock_header, "END");
        len = (hend-datablock_header)+80;

        // If BACKEND record is not present, insert it as first record.
        // TODO: Verify that we have room to insert the record.
        if(!ksearch(datablock_header, "BACKEND")) {
            // Move exsiting records to make room for new first record
            memmove(datablock_header+80, datablock_header, len);
            // Copy in BACKEND_RECORD string
            strncpy(datablock_header, BACKEND_RECORD, 80);
            // Increase len by 80 to account for the added record
            len += 80;
        }

        // Adjust length for any padding required for DirectIO
        if(directio) {
            // Round up to next multiple of 512
            len = (len+511) & ~511;
        }

        /* Write header (and padding, if any) */
        rv = write_all(fdraw, datablock_header, len);
        if (rv != len) {
            char msg[100];
            perror(thread_name);
            sprintf(msg, "Error writing data (datablock_header=%p, len=%d, rv=%d)", datablock_header, len, rv);
            hashpipe_error(thread_name, msg);
        }

        /* Write data */
        datablock_header = hpguppi_databuf_data(indb, curblock_in);
        len = blocksize;
        if(directio) {
            // Round up to next multiple of 512
            len = (len+511) & ~511;
        }
        rv = write_all(fdraw, datablock_header, (size_t)len);
        if (rv != len) {
            char msg[100];
            perror(thread_name);
            sprintf(msg, "Error writing data (datablock_header=%p, len=%d, rv=%d)", datablock_header, len, rv);
            hashpipe_error(thread_name, msg);
        }

        if(!directio) {
          /* flush output */
          fsync(fdraw);
        }

        /* Increment counter */
        block_count++;
      }
      /*** RAW Disk write out END*/
    }

    hpguppi_input_databuf_set_free(indb, curblock_in);
    blocks_per_second ++;

    // Update moving sum (for moving average)
    clock_gettime(CLOCK_MONOTONIC, &ts_free_input);
    fill_to_free_elapsed_ns = ELAPSED_NS(ts_block_recvd, ts_free_input);
    // Add new value, subtract old value
    fill_to_free_moving_sum_ns +=
        fill_to_free_elapsed_ns - fill_to_free_block_ns[curblock_in];
    // Store new value
    fill_to_free_block_ns[curblock_in] = fill_to_free_elapsed_ns;
    
    curblock_in  = (curblock_in + 1) % indb->header.n_block;

    /* Will exit if thread has been cancelled */
    pthread_testcancel();
  }
    
  hashpipe_info(thread_name, "exiting!");
  pthread_exit(NULL);

  pthread_cleanup_pop(0); /* Closes safe_close */
  pthread_cleanup_pop(0); /* Closes hpguppi_free_psrfits */
}

static hashpipe_thread_desc_t obs_rawdisk_thread = {
    name: "hpguppi_atasnap_obs_rawdisk_thread",
    skey: "OBSSTAT",
    init: NULL,
    run:  run,
    ibuf_desc: {hpguppi_input_databuf_create},
    obuf_desc: {NULL}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&obs_rawdisk_thread);
}

// vi: set ts=8 sw=4 et :
