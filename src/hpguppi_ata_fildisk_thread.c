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

#include "hpguppi_blade_databuf.h"
#include "hpguppi_params.h"
#include "hpguppi_udp.h"
#include "hpguppi_time.h"
#include "hpguppi_atasnap.h"
#include "hpguppi_rawspec.h"

#include "ioprio.h"

#include "mysigproc/mysigproc_utils.h"

// 80 character string for the BACKEND header record.
static const char BACKEND_RECORD[] =
// 0000000000111111111122222222223333333333
// 0123456789012345678901234567890123456789
  "BACKEND = 'FILTERBANK'                  " \
  "                                        ";

#define ELAPSED_S(start,stop) \
  ((int64_t)stop.tv_sec-start.tv_sec)

#define ELAPSED_NS(start,stop) \
  (ELAPSED_S(start,stop)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

void hpguppi_fil_read_header_from_status(
  char* statusbuf,
  sigproc_header_t* filheader
){
  rawspec_raw_hdr_t raw_hdr;
  rawspec_raw_parse_header(statusbuf, &raw_hdr);

  // raw_hdr.blocsize;
  // raw_hdr.npol;
  // raw_hdr.obsfreq;
  // raw_hdr.obsbw;
  // raw_hdr.directio;
  // raw_hdr.pktidx;
  // raw_hdr.nants;
  // raw_hdr.float_data;
  // raw_hdr.telescop;

  snprintf(filheader->source_name, sizeof(filheader->source_name), "%s", raw_hdr.src_name);

  filheader->machine_id = 0;                      // Machine ID
  hgetr8(statusbuf, "AZ", &filheader->az_start);  // azimuth
  hgetr8(statusbuf, "EL", &filheader->za_start);  // zenith angle
  filheader->za_start = 90.0 - filheader->za_start;
  filheader->nifs = raw_hdr.npol;                // Number of IF (stokes I => IF = 1) (AKA NPOL)

  filheader->telescope_id = fb_telescope_id(raw_hdr.telescop); // Telescope ID
  filheader->tstart = raw_hdr.mjd;    // MJD
  filheader->src_raj = 0.0;//raw_hdr.ra;    // RA in weird sigproc format
  filheader->src_dej = 0.0;//raw_hdr.dec;   // declination in sigproc format
  filheader->nbits = raw_hdr.nbits;   // Bit size
  filheader->nbeams = raw_hdr.nbeams; // Nbeams
  filheader->ibeam = raw_hdr.beam_id; // Beam number

  int nsources = raw_hdr.nants;
  if(filheader->nbeams > 0) {
    hashpipe_info(__FUNCTION__, "NBEAMS %d > 0 overriding NANTS %d value.", filheader->nbeams, raw_hdr.nants);
    nsources = filheader->nbeams;
  }

  // Emulate https://github.com/UCBerkeleySETI/rawspec/blob/162752186c8406c3ead985c3fd48f4035371cd0f/rawspec.c#L805-L830
  {
    // raw_hdr.obsnchan is total for all nants
    filheader->foff =
      raw_hdr.obsbw/(raw_hdr.obsnchan/nsources);
    // This computes first channel frequency (fch1).
    // raw_hdr.obsbw is always for single antenna
    // raw_hdr.obsnchan is total for all nants/nbeams
    filheader->fch1 = raw_hdr.obsfreq
      - raw_hdr.obsbw*((raw_hdr.obsnchan/nsources)-1)
          /(2*raw_hdr.obsnchan/nsources);
    
    filheader->nchans = raw_hdr.obsnchan/nsources; // Number of channels in file (across antenna/beams).
    filheader->tsamp = raw_hdr.tbin; // Sampling time in seconds.
  }
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
  uint32_t blocksize_perbeam=0, timesamples_perblock;
  int got_packet_0=0;
  fil_t** filbanks = NULL;
  sigproc_header_t filbank_header;

  unsigned char base_filename_stem_start;

  struct mjd_t *mjd = malloc(sizeof(struct mjd_t));

  /* Misc counters, etc */
  int i;

  uint64_t obs_npacket_total=0, obs_ndrop_total=0;
  uint32_t block_npacket=0, block_ndrop=0;

  uint64_t obs_start_pktidx = 0, obs_stop_pktidx = 0;
  uint64_t block_start_pktidx = 0, block_stop_pktidx = 0;
  
  char waiting=-1, flag_state_update=0;
  enum run_states state = IDLE;

  /* Filepath variables */
  char datadir[256];
  char fname[256];
  char *last_slash;

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
  hput_obsdone(st, 1);

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
      rv=hpguppi_databuf_wait_filled(indb, curblock_in);
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
    if (block_start_pktidx != block_stop_pktidx){
      hashpipe_warn(thread_name, "Current block seems out of order: starts at %lu, last ended at %lu.", block_start_pktidx, block_stop_pktidx);
    }
    hgetu8(datablock_header, "BLKSTOP", &block_stop_pktidx);

    switch(state_from_block_start_stop(obs_start_pktidx, obs_stop_pktidx, block_start_pktidx, block_stop_pktidx)){
      case IDLE:// If should IDLE, 
        if(state != IDLE){
          if(state == RECORD){//and recording, finalise block
            // If file open, close it
            if(filbanks) {
              for(i = 0; i < filbank_header.nbeams; i++){
                // Now destroy the fil struct
                // de-allocates memory and closes file
                if (destroy_fil(filbanks[i])) {
                  fprintf(stderr, "Error in destroy_fil#%d\n", i);
                  pthread_exit(NULL);
                }
                filbanks[i] = NULL;
              }
              // Reset filbank, got_packet_0
              free(filbanks);
              filbanks = NULL;
              got_packet_0 = 0;

              // Print end of recording conditions
              hashpipe_info(thread_name, "recording stopped: "
                "obs_start %lu obs_stop %lu blk_start_pktidx %lu blk_stop_pktidx %lu",
                obs_start_pktidx, obs_stop_pktidx, block_start_pktidx, block_stop_pktidx);
            }
            hput_obsdone(st, 1);
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
          hput_obsdone(st, 0);
          hgeti4(datablock_header, "STTVALID", &gp.stt_valid);
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
      /*** FIL Disk write out BEGIN */

      // Wait for packet 0 before starting write
      // "packet 0" is the first packet/block of the new recording,
      // it is not necessarily pktidx == 0.
      if (got_packet_0==0 && gp.stt_valid==1) {
        got_packet_0 = 1;
        // got_packet_0 only = 0 when filbanks is NULL, so don't bother clearing it
        hpguppi_fil_read_header_from_status(datablock_header, &filbank_header);
        filbanks = (fil_t**) malloc(filbank_header.nbeams*sizeof(fil_t*));

        // populate pf.basefilename //TODO this method is overkill for that purpose
        hpguppi_read_obs_params(datablock_header, &gp, &pf);

        // Create the output directory if needed        
        strncpy(datadir, pf.basefilename, 256);
        last_slash = strrchr(datadir, '/');
        if (last_slash!=NULL && last_slash!=datadir) {
          *last_slash = '\0';
          printf("Using directory '%s' for output.\n", datadir);
          if(mkdir_p(datadir, 0755) == -1) {
            hashpipe_error(thread_name, "mkdir_p(%s)", datadir);
            break;
          }
        }
        
        // finds last '/'
        base_filename_stem_start = strlen(pf.basefilename);
        while(base_filename_stem_start > 0 && pf.basefilename[base_filename_stem_start-1] != '/'){
          base_filename_stem_start--;
        }
        hashpipe_status_lock_safe(st);
          hputs(st->buf, "OBSSTEM", pf.basefilename+base_filename_stem_start);
        hashpipe_status_unlock_safe(st);

        // Open filterbank files for each beam
        for(i = 0; i < filbank_header.nbeams; i++) {
          sprintf(fname, "%s-beam%04d.fil", pf.basefilename, i);
          hashpipe_info(thread_name, "Opening filterbank file '%s'", fname);
          filbanks[i] = create_fil(fname, &rv, SIGPROC_CREATE_FIL_WRITE);
          if(filbanks[i]->file == NULL){
            hashpipe_error(thread_name, "Failed to open file: %s", fname);
          }
          memcpy(filbanks[i]->header, &filbank_header, sizeof(sigproc_header_t));
          filbanks[i]->header->ibeam = i;
          write_header(*filbanks[i]);
        }        

        hgetu4(datablock_header, "BLOCSIZE", &blocksize_perbeam);
        hgetu4(datablock_header, "PIPERBLK", &timesamples_perblock);
        blocksize_perbeam = blocksize_perbeam/filbank_header.nbeams;
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

        /* Write data */
        datablock_header = hpguppi_databuf_data(indb, curblock_in);
        for(i = 0; i < filbank_header.nbeams; i++) {
          rv = fwrite(
            datablock_header + i*blocksize_perbeam,
            1,
            blocksize_perbeam,
            filbanks[i]->file
          );
          if(rv != blocksize_perbeam
          ) {
              perror(thread_name);
              hashpipe_error(
                thread_name,
                "Error writing data (block#=%d, beam#=%d, len=%d, rv=%d)",
                curblock_in,
                i,
                blocksize_perbeam,
                rv
              );
              pthread_exit(NULL);
          }
        }
      }
      /*** FIL Disk write out END*/
    }

    hpguppi_databuf_set_free(indb, curblock_in);
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

static hashpipe_thread_desc_t fildisk_thread = {
    name: "hpguppi_ata_fildisk_thread",
    skey: "OBSSTAT",
    init: NULL,
    run:  run,
    ibuf_desc: {hpguppi_blade_output_databuf_create},
    obuf_desc: {NULL}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&fildisk_thread);
}

// vi: set ts=8 sw=4 et :
