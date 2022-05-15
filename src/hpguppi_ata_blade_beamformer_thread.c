#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <math.h>

#include <complex.h>

#include "hashpipe.h"
#include "hpguppi_time.h"
#include "hpguppi_databuf.h"
#include "hpguppi_blade_databuf.h"
#include "hpguppi_blade_ata_mode_b_capi.h"

#include "uvh5/uvh5_toml.h"
#include "radiointerferometryc99.h"
#include "antenna_weights.h"

#define ELAPSED_S(start,stop) \
  ((int64_t)stop.tv_sec-start.tv_sec)

#define ELAPSED_NS(start,stop) \
  (ELAPSED_S(start,stop)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

double mjd_mid_block(char* databuf_header){
  uint64_t pktidx, piperblk;
  struct mjd_t mjd = {0};
  uint64_t synctime = 0;
  double chan_bw = 1.0;

  double realtime_secs = 0.0;
  struct timespec ts;

  
  hgetu8(databuf_header, "PKTIDX", &pktidx);
  hgetu8(databuf_header, "PIPERBLK", &piperblk);
  hgetr8(databuf_header, "CHAN_BW", &chan_bw);
  hgetu8(databuf_header, "SYNCTIME", &synctime);

  // Calc real-time seconds since SYNCTIME for pktidx, taken to be a multiple of PKTNTIME:
  //
  //                          pktidx
  //     realtime_secs = -------------------
  //                        1e6 * chan_bw
  if(chan_bw != 0.0) {
    realtime_secs = (pktidx + piperblk/2) / (1e6 * fabs(chan_bw));
  }

  ts.tv_sec = (time_t)(synctime + rint(realtime_secs));
  ts.tv_nsec = (long)((realtime_secs - rint(realtime_secs)) * 1e9);

  get_mjd_from_timespec(&ts, &(mjd.stt_imjd), &(mjd.stt_smjd), &(mjd.stt_offs));
  //double t = ((double)mjd.stt_imjd) + ((double) mjd.stt_smjd + mjd.stt_offs)/(86400.0);
  //double t2 = ((double)synctime + (double)realtime_secs) / (double) 86400.0 + (double) 40587.0;
  //fprintf(stderr, "diff is: %.16f %.16f\n", t, t2);
  return (((double)mjd.stt_imjd) + ((double) mjd.stt_smjd + mjd.stt_offs)/(86400.0)) + 2400000.5;
}

// Populates `beam_coordinates` which should be nbeams*2 long (RA/DEC_OFFX pairs)
void collect_beamCoordinates(int nbeams, double* beam_coordinates, char* databuf_header) {
  char coordkey[9] = "DEC_OFFX";
  for(int beam_idx = 0; beam_idx < nbeams; beam_idx++) {
    sprintf(coordkey, "RA_OFF%d", beam_idx%10);
    hgetr8(databuf_header, coordkey, beam_coordinates+beam_idx*2+0);
    beam_coordinates[beam_idx*2+0] = calc_deg2rad(beam_coordinates[beam_idx*2+0]);

    sprintf(coordkey, "DEC_OFF%d", beam_idx%10);
    hgetr8(databuf_header, coordkey, beam_coordinates+beam_idx*2+1);
    beam_coordinates[beam_idx*2+1] = calc_deg2rad(beam_coordinates[beam_idx*2+1]);
  }
}


static void *run(hashpipe_thread_args_t *args)
{
  hpguppi_input_databuf_t *indb  = (hpguppi_input_databuf_t *)args->ibuf;
  hpguppi_blade_output_databuf_t *outdb = (hpguppi_blade_output_databuf_t *)args->obuf;
  char * databuf_header;

  hashpipe_status_t st = args->st;
  const char* status_key = args->thread_desc->skey;
  const char* thread_name = args->thread_desc->name;
  char buf_status[80];
  
  int curblock_in=0;
  int curblock_out=0;
  size_t i;

  /* Sundry flags */
  int hpguppi_databuf_wait_rv, status_state=0, update_status=1;

  /* Timestamp variables */
  struct timespec ts_status_update = {0}, ts_now = {0};
  const uint64_t status_update_period_ns = 1e9;
  
  uint64_t fill_to_free_elapsed_ns;
  uint64_t fill_to_free_moving_sum_ns = 0;
  uint64_t fill_to_free_block_ns[N_INPUT_BLOCKS] = {0};
  struct timespec ts_free_input = {0};
  struct timespec ts_blocks_recvd[N_INPUT_BLOCKS] = {0};

  char indb_data_dims_good_flag = 0;

  /* UVH5 header variables*/
	UVH5_header_t uvh5_header = {0};
  char tel_info_toml_filepath[70] = {'\0'};
  char obs_info_toml_filepath[70] = {'\0'};
  double obs_antenna_positions[BLADE_ATA_MODE_B_INPUT_NANT*3] = {0}, obs_beam_coordinates[BLADE_ATA_MODE_B_OUTPUT_NBEAM*2] = {0};
  struct blade_ata_mode_b_observation_meta observationMetaData = {0};
  struct LonLatAlt arrayReferencePosition = {0};
  
  /* BLADE variables */
  const size_t batch_size = 2;
  int input_output_blockid_pairs[N_INPUT_BLOCKS]; // input_id indexed associated output_id
  memset(input_output_blockid_pairs, -1, sizeof(int)*N_INPUT_BLOCKS);
  size_t dequeued_input_id = 0;

  double _Complex* antenna_calibration_coeffs;
  char obs_antenna_calibration_filepath[70] = {'\0'};
  char** obs_antenna_names = NULL;
  
  int cudaDeviceId = args->instance_id;
  
  hashpipe_status_lock_safe(&st);
  {
    hgeti4(st.buf, "CUDADEV", &cudaDeviceId);
  }
  hashpipe_status_unlock_safe(&st);
  if(cudaDeviceId >= 0){
    if(blade_use_device(cudaDeviceId)){
      hashpipe_info(args->thread_desc->name, "Successfully set CUDA device to %d.", cudaDeviceId);
    }
    else{
      hashpipe_info(args->thread_desc->name, "Failed to set CUDA device to %d.", cudaDeviceId);
      cudaDeviceId = -1;
    }
  }
  hashpipe_status_lock_safe(&st);
  {
    hputi4(st.buf, "CUDADEV", cudaDeviceId);
    hputi4(st.buf, "BLDBLKSZ", BLADE_BLOCK_DATA_SIZE);
  }
  hashpipe_status_unlock_safe(&st);

  //XXX Replace by actual calibration
  size_t size_of_calib = 20*192*2; //20 ants, 192chans, 2pol, HARCODED!!!!
  antenna_calibration_coeffs = 
    (double _Complex*) malloc(size_of_calib*sizeof (*antenna_calibration_coeffs));

  for (size_t i; i<size_of_calib; i++)
    antenna_calibration_coeffs[i] = 1.0 + 0.0*I;

  blade_ata_b_initialize(
    BLADE_ATA_MODE_B_CONFIG,
    batch_size,
    &observationMetaData,
    &arrayReferencePosition,
    obs_beam_coordinates,
    obs_antenna_positions,
    antenna_calibration_coeffs
  );

  if(BLADE_BLOCK_DATA_SIZE != blade_ata_b_get_output_size()*BLADE_ATA_MODE_B_OUTPUT_NCOMPLEX_BYTES){
    hashpipe_error(thread_name, "BLADE_BLOCK_DATA_SIZE %lu != %lu BLADE configured output size.", BLADE_BLOCK_DATA_SIZE, blade_ata_b_get_output_size()*BLADE_ATA_MODE_B_OUTPUT_NCOMPLEX_BYTES);
    pthread_exit(NULL);
  }

  for(i = 0; i < N_INPUT_BLOCKS; i++)
  {
    blade_pin_memory(hpguppi_databuf_data(indb, i), BLOCK_DATA_SIZE);
    blade_pin_memory(hpguppi_databuf_data(outdb, i), BLADE_BLOCK_DATA_SIZE);
  }

  int32_t input_buffer_dim_NANTS = 0, prev_flagged_NANTS = 0;
  int32_t input_buffer_dim_NCHAN = 0, prev_flagged_NCHAN = 0;
  int32_t input_buffer_dim_NTIME = 0, prev_flagged_NTIME = 0;
  int32_t input_buffer_dim_NPOLS = 0, prev_flagged_NPOLS = 0;

  int64_t pktidx_obs_start, pktidx_obs_start_prev, pktidx_blk_start;
  int fenchan;
  double timemjd_midblock=0, dut1=0; // DUT1
  timemjd_midblock += 1;
  double timejd_midblock=0;

  while (run_threads())
  {
    do {
      hpguppi_databuf_wait_rv = hpguppi_databuf_wait_filled(
          indb, curblock_in);
      
      clock_gettime(CLOCK_MONOTONIC, &ts_now);

      update_status = ELAPSED_NS(ts_status_update, ts_now) > status_update_period_ns;

      if(hpguppi_databuf_wait_rv == HASHPIPE_TIMEOUT && !update_status) {
        // No, continue receiving
        continue;
      }

      // We perform some status buffer updates every second
      if(update_status) {
        sprintf(buf_status, "%d/%d", hpguppi_databuf_total_status(outdb), N_INPUT_BLOCKS);

        hashpipe_status_lock_safe(&st);
        {
          hputi4(st.buf, "BLDBLKSZ", BLADE_BLOCK_DATA_SIZE);
          hputs(st.buf, "BLDBUFST", buf_status);
          hputr4(st.buf, "BLDBLKMS",
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
    
    indb_data_dims_good_flag = 1; // innocent until proven guilty
    databuf_header = hpguppi_databuf_header(indb, curblock_in);
    hgeti4(databuf_header, "NANTS", &input_buffer_dim_NANTS);
    hgeti4(databuf_header, "NCHAN", &input_buffer_dim_NCHAN);
    hgeti4(databuf_header, "PIPERBLK", &input_buffer_dim_NTIME);
    hgeti4(databuf_header, "NPOL", &input_buffer_dim_NPOLS);

    if(input_buffer_dim_NANTS != BLADE_ATA_MODE_B_INPUT_NANT){
      indb_data_dims_good_flag = 0;
      if(prev_flagged_NANTS != input_buffer_dim_NANTS){
        prev_flagged_NANTS = input_buffer_dim_NANTS;
        hashpipe_error(thread_name, "Incoming data_buffer has NANTS %lu != %lu. Ignored.", input_buffer_dim_NANTS, BLADE_ATA_MODE_B_INPUT_NANT);
      }
    }
    else if(input_buffer_dim_NCHAN != BLADE_ATA_MODE_B_ANT_NCHAN){
      indb_data_dims_good_flag = 0;
      if(prev_flagged_NCHAN != input_buffer_dim_NCHAN){
        prev_flagged_NCHAN = input_buffer_dim_NCHAN;
        hashpipe_error(thread_name, "Incoming data_buffer has NCHANS %lu != %lu. Ignored.", input_buffer_dim_NCHAN, BLADE_ATA_MODE_B_ANT_NCHAN);
      }
    }
    else if(input_buffer_dim_NTIME != BLADE_ATA_MODE_B_NTIME){
      indb_data_dims_good_flag = 0;
      if(prev_flagged_NTIME != input_buffer_dim_NTIME){
        prev_flagged_NTIME = input_buffer_dim_NTIME;
        hashpipe_error(thread_name, "Incoming data_buffer has NTIME %lu != %lu. Ignored.", input_buffer_dim_NTIME, BLADE_ATA_MODE_B_NTIME);
      }
    }
    else if(input_buffer_dim_NPOLS != BLADE_ATA_MODE_B_NPOL){
      indb_data_dims_good_flag = 0;
      if(prev_flagged_NPOLS != input_buffer_dim_NPOLS){
        prev_flagged_NPOLS = input_buffer_dim_NPOLS;
        hashpipe_error(thread_name, "Incoming data_buffer has NPOLS %lu != %lu. Ignored.", input_buffer_dim_NPOLS, BLADE_ATA_MODE_B_NPOL);
      }
    }

    if (!indb_data_dims_good_flag) {
      hpguppi_databuf_set_free(indb, curblock_in);
      curblock_in  = (curblock_in + 1) % indb->header.n_block;
      continue;
    }
    else{
      prev_flagged_NANTS = input_buffer_dim_NANTS;
      prev_flagged_NCHAN = input_buffer_dim_NCHAN;
      prev_flagged_NTIME = input_buffer_dim_NTIME;
      prev_flagged_NPOLS = input_buffer_dim_NPOLS;
    }

    hgeti8(databuf_header, "PKTSTART", &pktidx_obs_start);
    hgeti8(databuf_header, "BLKSTART", &pktidx_blk_start);

    if(pktidx_obs_start != pktidx_obs_start_prev && pktidx_obs_start >= pktidx_blk_start){ // first block of observation
        hgetu8(databuf_header, "SCHAN", &observationMetaData.frequencyStartIndex);
        hgetr8(databuf_header, "CHAN_BW", &observationMetaData.channelBandwidthHz);
        hgeti4(databuf_header, "FENCHAN", &fenchan);
        hgetr8(databuf_header, "OBSFREQ", &observationMetaData.rfFrequencyHz);

        double tmp = (double)observationMetaData.rfFrequencyHz +
          (-(double)observationMetaData.frequencyStartIndex - ((double)input_buffer_dim_NCHAN / 2.0)
            + ((double)fenchan / 2.0) + 0.5
          ) * (double)observationMetaData.channelBandwidthHz;


        observationMetaData.rfFrequencyHz = tmp;

        observationMetaData.rfFrequencyHz *= 1e6;
        observationMetaData.channelBandwidthHz *= 1e6;
        observationMetaData.totalBandwidthHz = fenchan * observationMetaData.channelBandwidthHz;

        if(obs_antenna_names){
          free(obs_antenna_names);
        }

        hashpipe_status_lock_safe(&st);
        {
          hgets(st.buf, "TELINFOP", 70, tel_info_toml_filepath);
          hgets(st.buf, "OBSINFOP", 70, obs_info_toml_filepath);
          hgets(st.buf, "CALWGHTP", 70, obs_antenna_calibration_filepath);
        }
        hashpipe_status_unlock_safe(&st);

        hashpipe_info(thread_name, "Parsing '%s' as Telescope information.", tel_info_toml_filepath);
        UVH5toml_parse_telescope_info(tel_info_toml_filepath, &uvh5_header);
        hashpipe_info(thread_name, "Parsing '%s' as Observation information.", obs_info_toml_filepath);
        UVH5toml_parse_observation_info(obs_info_toml_filepath, &uvh5_header);
        UVH5Hadmin(&uvh5_header);
        arrayReferencePosition.LAT = calc_deg2rad(uvh5_header.latitude);
        arrayReferencePosition.LON = calc_deg2rad(uvh5_header.longitude);
        arrayReferencePosition.ALT = uvh5_header.altitude;

        obs_antenna_names = malloc(uvh5_header.Nants_data*sizeof(char*));

        // At this point we have XYZ uvh5_header.antenna_positions, and ENU uvh5_header._antenna_enu_positions
        for(i = 0; i < uvh5_header.Nants_data; i++){
          int ant_idx = uvh5_header._antenna_num_idx_map[
            uvh5_header._antenna_numbers_data[i]
          ];
          obs_antenna_positions[i*3 + 0] = uvh5_header.antenna_positions[ant_idx*3 + 0];
          obs_antenna_positions[i*3 + 1] = uvh5_header.antenna_positions[ant_idx*3 + 1];
          obs_antenna_positions[i*3 + 2] = uvh5_header.antenna_positions[ant_idx*3 + 2];

          obs_antenna_names[i] = uvh5_header.antenna_names[ant_idx];
        }
        // BLADE needs ECEF coordinates
        // It doesn't make sense to convert from ECEF back to XYZ in the uvh5 library
        // but then back to ECEF here. TODO don't do the above
        calc_position_to_ecef_frame_from_xyz(
          obs_antenna_positions,
          uvh5_header.Nants_data,
          arrayReferencePosition.LON,
          arrayReferencePosition.LAT,
          arrayReferencePosition.ALT
        );
        // observationMetaData.referenceAntennaIndex = uvh5_header._antenna_num_idx_map[
        //   uvh5_header._antenna_numbers_data[0]
        // ];
        observationMetaData.referenceAntennaIndex = 0;

        free(antenna_calibration_coeffs);
        collect_beamCoordinates(BLADE_ATA_MODE_B_OUTPUT_NBEAM,
            obs_beam_coordinates, obs_phase_center, databuf_header);
        if(
          read_antenna_weights(
            obs_antenna_calibration_filepath,
            uvh5_header.Nants_data, // number of antenna of interest
            obs_antenna_names, // the antenna of interest
            observationMetaData.frequencyStartIndex, // the first channel
            input_buffer_dim_NCHAN, // the number of channels
            &antenna_calibration_coeffs // return value
          )
        ){
          // Failed to open CALWGHTP file, set 1.0+0.0j
          hashpipe_warn(thread_name, "CALWGHTP `%s` could not be opened. Using 1.0+0.0j.", obs_antenna_calibration_filepath);
          errno = 0;
          size_of_calib = 
            input_buffer_dim_NANTS*
            input_buffer_dim_NCHAN*
            input_buffer_dim_NPOLS;
          antenna_calibration_coeffs = malloc(size_of_calib*sizeof(double _Complex*));

          for(i = 0; i < size_of_calib; i++){
            antenna_calibration_coeffs[i] = 1.0 + 0.0*I;
          }
        }

        // Free previously queued buffers
        for(i = 0; i < indb->header.n_block; i++)
        {
          if(input_output_blockid_pairs[i] != -1) {
            hpguppi_databuf_set_free(indb, i);
            // set filled old output, to synchronise the disk thread
            hpguppi_databuf_set_filled(outdb, input_output_blockid_pairs[i]);
            input_output_blockid_pairs[i] = -1;
          }
        }

        blade_ata_b_terminate();
        blade_ata_b_initialize(
          BLADE_ATA_MODE_B_CONFIG,
          batch_size,
          &observationMetaData,
          &arrayReferencePosition,
          obs_beam_coordinates,
          obs_antenna_positions,
          antenna_calibration_coeffs
        );
    }
    pktidx_obs_start_prev = pktidx_obs_start;
    hgetr8(databuf_header, "DUT1", &dut1);
    timemjd_midblock = mjd_mid_block(databuf_header);


    // waiting for output buffer to be free
    while ((hpguppi_databuf_wait_rv=hpguppi_databuf_wait_free(outdb, curblock_out)) !=
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
        hashpipe_error(thread_name, "error waiting for output buffer #%d, rv: %i", curblock_out, hpguppi_databuf_wait_rv);
        pthread_exit(NULL);
        break;
      }
    }
    // collect free output_buffer 
    if(hpguppi_databuf_wait_rv != HASHPIPE_OK){
      // presume an error occurred while status_state for a free output buffer
      break;
    }

    if(status_state != 3){
      hashpipe_status_lock_safe(&st);
      hputs(st.buf, status_key, "beamforming");
      hashpipe_status_unlock_safe(&st);
      status_state = 3;
    }

    // Asynchronously queue
    input_output_blockid_pairs[curblock_in] = curblock_out;
    while(!
      blade_ata_b_enqueue(
        (void*) hpguppi_databuf_data(indb, curblock_in),
        (void*) hpguppi_databuf_data(outdb, curblock_out),
        curblock_in,
        timejd_midblock,
        dut1
      )
    ){
      if(status_state != 4){
        hashpipe_status_lock_safe(&st);
        hputs(st.buf, status_key, "beamforming blocked");
        hashpipe_status_unlock_safe(&st);
        status_state = 4;
      } 
    }
    clock_gettime(CLOCK_MONOTONIC, &ts_blocks_recvd[curblock_in]);

    {// Asynchronous CPU work 
      // copy across the header
      databuf_header = hpguppi_databuf_header(outdb, curblock_out);
      memcpy(databuf_header, 
            hpguppi_databuf_header(indb, curblock_in), 
            BLOCK_HDR_SIZE);
      
      //TODO upate output_buffer headers to reflect that they contain beams
      hputi4(databuf_header, "OBSNCHAN", BLADE_ATA_MODE_B_ANT_NCHAN); // beams are split into separate files...
      hputi4(databuf_header, "NBEAMS", BLADE_ATA_MODE_B_OUTPUT_NBEAM);
      hputi4(databuf_header, "NBITS", BLADE_ATA_MODE_B_OUTPUT_NCOMPLEX_BYTES*8/2);
      hputs(databuf_header, "DATATYPE", "FLOAT");
      hputs(databuf_header, "SMPLTYPE", BLADE_ATA_MODE_B_OUTPUT_NCOMPLEX_BYTES == 8 ? "CF32" : "CF16");
      hputi4(databuf_header, "BLOCSIZE", BLADE_BLOCK_DATA_SIZE);
    }

    // hashpipe_info(thread_name, "batched block #%d as input #%d.", curblock_in, inputs_batched);
    curblock_in  = (curblock_in + 1) % indb->header.n_block;
    // hashpipe_info(thread_name, "batched block #%d as output #%d.", curblock_out, outputs_batched);
    curblock_out = (curblock_out + 1) % outdb->header.n_block;

    // Dequeue all completed buffers
    while (blade_ata_b_dequeue(&dequeued_input_id))
    {
      hpguppi_databuf_set_free(indb, dequeued_input_id);
      hpguppi_databuf_set_filled(outdb, input_output_blockid_pairs[dequeued_input_id]);
      input_output_blockid_pairs[dequeued_input_id] = -1;

      // Update moving sum (for moving average)
      clock_gettime(CLOCK_MONOTONIC, &ts_free_input);
      fill_to_free_elapsed_ns = ELAPSED_NS(ts_blocks_recvd[dequeued_input_id], ts_free_input)/batch_size;
      fill_to_free_moving_sum_ns +=
          fill_to_free_elapsed_ns - fill_to_free_block_ns[dequeued_input_id];
      // Store new value
      fill_to_free_block_ns[dequeued_input_id] = fill_to_free_elapsed_ns;
    }
  }

  hashpipe_info(thread_name, "returning");
  blade_ata_b_terminate();
  return NULL;
}

static hashpipe_thread_desc_t blade_beamformer_thread = {
  name: "hpguppi_ata_blade_beamformer_thread",
  skey: "BEAMSTAT",
  init: NULL,
  run: run,
  ibuf_desc: {hpguppi_input_databuf_create},
  obuf_desc: {hpguppi_blade_output_databuf_create}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(& blade_beamformer_thread);
}
