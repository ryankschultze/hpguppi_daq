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

#include "hpguppi_xgpu_databuf.h"
#include "hpguppi_params.h"
#include "hpguppi_udp.h"
#include "hpguppi_time.h"
#include "hpguppi_atasnap.h"

#include "ioprio.h"

#include "uvh5.h"
#include "uvh5/uvh5_toml.h"
#include "uvh5/uvh5_bool_t.h"
#include "radiointerferometryc99.h"

#define MJD0 2400000.5

#define ELAPSED_S(start,stop) \
  ((int64_t)stop.tv_sec-start.tv_sec)

#define ELAPSED_NS(start,stop) \
  (ELAPSED_S(start,stop)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

static int safe_close(UVH5_file_t* uvh5_file) {
  UVH5close(uvh5_file);
  return 0;
}

#define HPUT_DAQ_STATE(st, state)\
  hputs(st->buf, "DAQSTATE", state == IDLE  ? "idling" :\
                             state == ARMED ? "armed"  :\
                             "recording")

static void *run(hashpipe_thread_args_t * args)
{
    // Local aliases to shorten access to args fields
    // Our output buffer happens to be a hpguppi_xgpu_output_databuf

  hpguppi_output_xgpu_databuf_t* indb = (hpguppi_output_xgpu_databuf_t *) args->ibuf;

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

	UVH5_file_t uvh5_file = {0};
	UVH5_header_t* uvh5_header = &uvh5_file.header;
  pthread_cleanup_push((void *)safe_close, &uvh5_file);

  /* Set I/O priority class for this thread to "real time" */
  if(ioprio_set(IOPRIO_WHO_PROCESS, 0, IOPRIO_PRIO_VALUE(IOPRIO_CLASS_RT, 7))) {
    hashpipe_error(thread_name, "ioprio_set IOPRIO_CLASS_RT");
  }

  char history[] =
    "hpguppi_daq: plug_here\n"
    "uvh5c99: https://github.com/MydonSolutions/uvh5c99\n";

  int rv = 0;
  int curblock_in=0;

  char *datablock_header;
  int got_block_0=0;

  unsigned char base_filename_stem_start;

  struct mjd_t *mjd = malloc(sizeof(struct mjd_t));
  double longitude_rad, latitude_rad, hour_angle_rad, declination_rad;
  double ra_rad, dec_rad;

  double dut1 = 0.0;
  double tau;
  double chan_bw, obs_freq;
  size_t xgpu_output_elements;

  /* Misc counters, etc */
  int i, j;
  uint64_t obs_npacket_total=0, obs_ndrop_total=0;
  uint64_t ndrop_obs_start=0, ndrop_obs_current=0;
  uint32_t block_npacket=0, block_ndrop=0;

  uint64_t obs_start_pktidx = 0, obs_stop_pktidx = 0;
  uint64_t block_start_pktidx = 0, block_stop_pktidx = 0;

  char waiting=-1, flag_state_update=0;
  enum run_states state = IDLE;

  int32_t stt_valid;

  // Used to calculate moving average of fill-to-free times for input blocks
  uint64_t fill_to_free_elapsed_ns;
  uint64_t fill_to_free_moving_sum_ns = 0;
  uint64_t fill_to_free_block_ns[N_XGPU_OUTPUT_BLOCKS] = {0};
  struct timespec ts_free_input = {0}, ts_block_recvd = {0};

  /* Heartbeat variables */
  time_t lasttime = 0;
  time_t curtime = 0;
  char timestr[32] = {0};

  char tel_info_toml_filepath[70] = {'\0'};
  char obs_info_toml_filepath[70] = {'\0'};

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
              hgetu8(st->buf, "NDROP", &ndrop_obs_current);
              hputu8(st->buf, "OBSNPKTS", obs_npacket_total);
              hputu8(st->buf, "OBSNDROP", ndrop_obs_current - ndrop_obs_start);
              hputu4(st->buf, "OBSBLKPS", blocks_per_second);
              hputr4(st->buf, "OBSBLKMS",
                round((double)fill_to_free_moving_sum_ns / N_XGPU_OUTPUT_BLOCKS) / 1e6);
              hputs(st->buf, "DAQPULSE", timestr);
              HPUT_DAQ_STATE(st, state);
          }
          hashpipe_status_unlock_safe(st);
          blocks_per_second = 0;
      }

      // Waiting for input
      rv=hpguppi_output_xgpu_databuf_wait_filled(indb, curblock_in);
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

    datablock_header = hpguppi_xgpu_output_databuf_header(indb, curblock_in);
    hgetu8(datablock_header, "PKTSTART", &obs_start_pktidx);
    hgetu8(datablock_header, "PKTSTOP", &obs_stop_pktidx);
    hgetu8(datablock_header, "BLKSTART", &block_start_pktidx);
    // if (block_start_pktidx != block_stop_pktidx){
    //   hashpipe_warn(thread_name, "Current block seems out of order: starts at %lu, last ended at %lu.", block_start_pktidx, block_stop_pktidx);
    // }
    hgetu8(datablock_header, "BLKSTOP", &block_stop_pktidx);

    switch(state_from_block_start_stop(obs_start_pktidx, obs_stop_pktidx, block_start_pktidx, block_stop_pktidx)){
      case IDLE:// If should IDLE,
        if(state != IDLE){
          if(state == RECORD){//and recording, finalise block
            // If file open, close it
            if(uvh5_file.file_id) {
              // Close file
              free(uvh5_header->object_name);
              UVH5close(&uvh5_file);
              memset(&uvh5_file, 0, sizeof(UVH5_file_t));
              uvh5_header = &uvh5_file.header;
              got_block_0 = 0;
            }
          }
          // Print end of recording conditions
          hashpipe_info(thread_name, "%s ended: "
            "obs_start %lu obs_stop %lu blk_start_pktidx %lu blk_stop_pktidx %lu",
            state == RECORD ? "recording" : "arming",
            obs_start_pktidx, obs_stop_pktidx, block_start_pktidx, block_stop_pktidx);
          hput_obsdone(st, 1);
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
            state = ARMED;
            update_stt_status_keys(st, state, obs_start_pktidx, mjd);
            hputu4(datablock_header, "STTVALID", 1);
            hputu4(datablock_header, "STT_IMJD", mjd->stt_imjd);
            hputu4(datablock_header, "STT_SMJD", mjd->stt_smjd);
            hputr8(datablock_header, "STT_OFFS", mjd->stt_offs);
            hput_obsdone(st, 0);
          }
          flag_state_update = 1;
          state = RECORD;
          hgeti4(datablock_header, "STTVALID", &stt_valid);
          
          hashpipe_status_lock_safe(st);
          {
            hgetu8(st->buf, "NDROP", &ndrop_obs_start);
          }
          hashpipe_status_unlock_safe(st);
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
          hput_obsdone(st, 0);
        }
      default:
        break;
    }

    if(state == RECORD){
      /***  Disk write out BEGIN */

      // Wait for block 0 before starting write
      // "block 0" is the first block of the new recording
      if (got_block_0==0 && stt_valid == 1) {
        got_block_0 = 1;

        { // Setup UVH5 File object
          // assume everything is closed, start up for new
          // set header scalar data
          hgeti4(datablock_header, "NCHAN", &uvh5_header->Nfreqs);
          hgetr8(datablock_header, "CHAN_BW", &chan_bw);
          hgetr8(datablock_header, "OBSFREQ", &obs_freq);
          hgets(datablock_header, "UVH5TELP", 70, tel_info_toml_filepath);
          hgets(datablock_header, "UVH5OBSP", 70, obs_info_toml_filepath);
          uvh5_header->Nspws = 1;
          uvh5_header->Ntimes = 0; // initially
          uvh5_header->Nblts = 0; // uvh5_header->Nbls * uvh5_header->Ntimes;

          hashpipe_info(thread_name, "Parsing '%s' as Telescope information.", tel_info_toml_filepath);
          UVH5toml_parse_telescope_info(tel_info_toml_filepath, uvh5_header);
          hashpipe_info(thread_name, "Parsing '%s' as Observation information.", obs_info_toml_filepath);
          UVH5toml_parse_observation_info(obs_info_toml_filepath, uvh5_header);
          UVH5Hadmin(uvh5_header);

          uvh5_header->spw_array[0] = 1;
          obs_freq -= uvh5_header->Nfreqs*chan_bw/2;
          uvh5_header->channel_width[0] = chan_bw * 1e6;
          for(i = 0; i < uvh5_header->Nfreqs; i++) {
            uvh5_header->channel_width[i] = uvh5_header->channel_width[0];
            uvh5_header->freq_array[i] = (obs_freq + (i+0.5)*chan_bw) * 1e6;
          }

          uvh5_header->instrument = "ATA";
          uvh5_header->object_name = malloc(71);
          uvh5_header->object_name[70] = '\0';
          hgets(datablock_header, "SRC_NAME", 70, uvh5_header->object_name);

          uvh5_header->history = history;
          uvh5_header->phase_type = "phased";
          hgetr8(datablock_header, "RA_STR", &uvh5_header->phase_center_ra);
          hgetr8(datablock_header, "DEC_STR", &uvh5_header->phase_center_dec);
          uvh5_header->phase_center_ra = calc_deg2rad(uvh5_header->phase_center_ra);
          uvh5_header->phase_center_dec = calc_deg2rad(uvh5_header->phase_center_dec);
          uvh5_header->phase_center_epoch = 2000.0;
          uvh5_header->phase_center_frame = "icrs";

          longitude_rad = calc_deg2rad(uvh5_header->longitude);
          latitude_rad = calc_deg2rad(uvh5_header->latitude);

          dut1 = 0.0;
          hgetr8(datablock_header, "DUT1", &dut1); // single DUT1 value for all observation time
          hgetr8(datablock_header, "XTIMEINT", &tau); // calculated in xgpu_thread from tbin*nsamperblk*blks_per_integration
          hashpipe_info(thread_name, "DUT1: %f", dut1);

          // uvh5_header->lst_array = malloc(sizeof(double) * uvh5_header->Nbls);

          uvh5_header->time_array[0] = MJD0;
          uvh5_header->time_array[0] += (double)mjd->stt_imjd;
          uvh5_header->time_array[0] += ((double)mjd->stt_smjd)/86400.0;
          uvh5_header->time_array[0] += (tau/2) / DAYSEC;
          hashpipe_info(thread_name, "UVH5 time starts at %f (tau/DAYSEC %f)", uvh5_header->time_array[0], tau/DAYSEC);
          uvh5_header->integration_time[0] = tau;
          // uvh5_header->lst_array[0] = calc_lst(uvh5_header->time_array[0], dut1) + longitude_rad;
          uvh5_header->dut1 = dut1;
          for (i = 1; i < uvh5_header->Nbls; i++)
          {
            uvh5_header->time_array[i] = uvh5_header->time_array[0];
            // uvh5_header->lst_array[i] = uvh5_header->lst_array[0];
            uvh5_header->integration_time[i] = uvh5_header->integration_time[0];
          }

          hgetu8(datablock_header, "X_TRILEN", &xgpu_output_elements);
          hashpipe_info(thread_name, 
            "Picked up xgpu output element count of %llu (%llu products per %d frequencies)", 
            xgpu_output_elements, xgpu_output_elements/uvh5_header->Nfreqs, uvh5_header->Nfreqs
          );
        }// Setup UVH5 File object
        
        hpguppi_read_obs_params(datablock_header, &gp, &pf);
        hpguppi_read_subint_params(datablock_header, &gp, &pf);
        char fname[256];
        sprintf(fname, "%s.uvh5", pf.basefilename);
        fprintf(stderr, "Opening uvh5 file '%s'\n", fname);

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
        UVH5open(fname, &uvh5_file, UVH5TcreateCI32()); // TODO if DP4A
        // if () {
        //   hashpipe_error(thread_name, "Error opening file.");
        //   pthread_exit(NULL);
        // }
      }

      /* If we got block 0, write data to disk */
      if (got_block_0) {
        if(waiting != -1){
          /* Note writing status */
          waiting = -1;
          hashpipe_status_lock_safe(st);
          hputs(st->buf, status_key, "writing");
          hashpipe_status_unlock_safe(st);
        }

        // memset(uvh5_file.visdata, 1, uvh5_header->Nbls*uvh5_header->Npols*uvh5_header->Nfreqs);
        UVH5visdata_from_xgpu_int_output(
          (UVH5_CI32_t*) hpguppi_xgpu_output_databuf_data(indb, curblock_in),
          (UVH5_CI32_t*) uvh5_file.visdata,
          xgpu_output_elements,
          &uvh5_file.header
        );

        uvh5_header->time_array[0] += + tau/DAYSEC;
        // uvh5_header->lst_array[0] = calc_lst(uvh5_header->time_array[0], dut1) + longitude_rad;
        hgetr8(datablock_header, "NSAMPLES", uvh5_file.nsamples);
        for (i = 0; i < uvh5_header->Nbls; i++) {
          uvh5_header->time_array[i] = uvh5_header->time_array[0];
          // uvh5_header->lst_array[i] = uvh5_header->lst_array[0];
          for (j = 0; j < uvh5_header->Nfreqs*uvh5_header->Npols; j++) {
            uvh5_file.nsamples[i*uvh5_header->Nfreqs*uvh5_header->Npols + j] = uvh5_file.nsamples[0];
            uvh5_file.flags[i*uvh5_header->Nfreqs*uvh5_header->Npols + j] = UVH5_FALSE;
          }
        }

        if (1) { // Phased
          hgetr8(datablock_header, "RA_STR", &ra_rad);
          hgetr8(datablock_header, "DEC_STR", &dec_rad);
          ra_rad = calc_deg2rad(ra_rad);
          dec_rad = calc_deg2rad(dec_rad);

          calc_ha_dec_rad(
            ra_rad,
            dec_rad,
            longitude_rad,
            latitude_rad,
            uvh5_header->altitude,
            uvh5_header->time_array[0],
            dut1,
            &hour_angle_rad,
            &declination_rad
          );

          memcpy(uvh5_header->_antenna_uvw_positions, uvh5_header->_antenna_enu_positions, sizeof(double)*uvh5_header->Nants_telescope*3);
          calc_position_to_uvw_frame_from_enu(
            uvh5_header->_antenna_uvw_positions,
            uvh5_header->Nants_telescope,
            hour_angle_rad,
            declination_rad,
            latitude_rad
          );

          UVH5permutate_uvws(uvh5_header);
        }

        UVH5write_dynamic(&uvh5_file);
      }

      // If obviously last block, close up and transition to IDLE
      if(block_stop_pktidx >= obs_stop_pktidx){
        // If file open, close it
        if(uvh5_file.file_id) {
          // Close file
          free(uvh5_header->object_name);
          UVH5close(&uvh5_file);
          memset(&uvh5_file, 0, sizeof(UVH5_file_t));
          uvh5_header = &uvh5_file.header;
          got_block_0 = 0;

          // Print end of recording conditions
          hashpipe_info(thread_name, "recorded last block: "
            "obs_start %lu obs_stop %lu blk_start_pktidx %lu blk_stop_pktidx %lu",
            obs_start_pktidx, obs_stop_pktidx, block_start_pktidx, block_stop_pktidx);
        }
        hput_obsdone(st, 1);
        flag_state_update = 1;
        state = IDLE;
        // Update STT with state = IDLE, setting STTVALID = 0
        update_stt_status_keys(st, state, obs_start_pktidx, mjd);
      }
      /*** UVH5 Disk write out END*/
    }

    hpguppi_output_xgpu_databuf_set_free(indb, curblock_in);
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

static hashpipe_thread_desc_t obs_uvh5disk_thread = {
    name: "hpguppi_ata_obs_uvh5disk_thread",
    skey: "OBSSTAT",
    init: NULL,
    run:  run,
    ibuf_desc: {hpguppi_output_xgpu_databuf_create},
    obuf_desc: {NULL}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&obs_uvh5disk_thread);
}

// vi: set ts=8 sw=4 et :
