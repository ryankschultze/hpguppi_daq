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

#define HPGUPPI_DAQ_CONTROL "/tmp/hpguppi_daq_control"
#define MAX_CMD_LEN 1024

#define PKTSOCK_BYTES_PER_FRAME (16384)
#define PKTSOCK_FRAMES_PER_BLOCK (8)
#define PKTSOCK_NBLOCKS (800)
#define PKTSOCK_NFRAMES (PKTSOCK_FRAMES_PER_BLOCK * PKTSOCK_NBLOCKS)

#define HPUT_DAQ_STATE(st, state)\
  hputs(st->buf, "DAQSTATE", state == IDLE  ? "idling" :\
                             state == ARMED ? "armed"  :\
                             "recording")
// Define run states.  Currently three run states are defined: IDLE, ARMED, and
// RECORD.
//
// In all states, the PKTIDX field is updated with the value from received
// packets.  Whenever the first PKTIDX of a block is received (i.e. whenever
// PKTIDX is a multiple of the number of packets per block), the value for
// PKTSTART and DWELL are read from the status buffer.  PKTSTART is rounded
// down, if needed, to ensure that it is a multiple of the number of packets
// per block, then PKTSTART is written back to the status buffer.  DWELL is
// interpreted as the number of seconds to record and is used to calculate
// PKTSTOP (which gets rounded down, if needed, to be a multiple of the number
// of packets per block).
//
// In the IDLE state, incoming packets are dropped, but PKTIDX is still
// updated.  Transitioning out of the IDLE state will reinitialize the current
// blocks.  A "start" command will transition the state from IDLE to ARMED.
// A "stop" command will transition the state from ARMED or RECORD to IDLE.
//
// In the ARMED state, incoming packets are processed (i.e. stored in the net
// thread's output buffer) and full blocks are passed to the next thread.  When
// the processed PKTIDX is equal to PKTSTART the state transitions to RECORD
// and the following actions occur:
//
//   1. The MJD of the observation start time is calculated (TODO Calculate
//      from SYNCTIME and PKTIDX)
//
//   2. The packet stats counters are reset
//
//   3. The STT_IMDJ and STT_SMJD are updated in the status buffer
//
//   4. STTVALID is set to 1
//
// In the RECORD state, incoming packets are processed (i.e. stored in the net
// thread's output buffer) and full blocks are passed to the next thread (same
// as in the ARMED state).  When the processed PKTIDX is greater than or equal
// to PKTSTOP the state transitions to ARMED and STTVALID is set to 0.
//
// The downstream thread (i.e. hpguppi_rawdisk_thread) is expected to use a
// combination of PKTIDX, PKTSTART, PKTSTOP, and (optionally) STTVALID to
// determine whether the blocks should be discarded or processed (e.g. written
// to disk).

enum run_states {IDLE, ARMED, RECORD};

//  if state == RECORD && STTVALID == 0 
//    STTVALID=1
//    calculate and store STT_IMJD, STT_SMJD
//  else if STTVALID != 0
//    STTVALID=0
//  endif
static void update_stt_status_keys( hashpipe_status_t *st,
                                    enum run_states state,
                                    uint64_t pktidx,
                                    struct mjd_t *mjd){
  // uint32_t pktntime = ATASNAP_DEFAULT_PKTNTIME;
  uint64_t synctime = 0;
  double chan_bw = 1.0;

  double realtime_secs = 0.0;
  struct timespec ts;

  uint32_t sttvalid = 0;
  hashpipe_status_lock_safe(st);
  {
    hgetu4(st->buf, "STTVALID", &sttvalid);
    if((state == ARMED || state == RECORD) && sttvalid != 1) {
      hputu4(st->buf, "STTVALID", 1);

      // hgetu4(st->buf, "PKTNTIME", &pktntime);
      hgetr8(st->buf, "CHAN_BW", &chan_bw);
      hgetu8(st->buf, "SYNCTIME", &synctime);

      // Calc real-time seconds since SYNCTIME for pktidx, taken to be a multiple of PKTNTIME:
      //
      //                          pktidx
      //     realtime_secs = -------------------
      //                        1e6 * chan_bw
      if(chan_bw != 0.0) {
        realtime_secs = pktidx / (1e6 * fabs(chan_bw));
      }

      ts.tv_sec = (time_t)(synctime + rint(realtime_secs));
      ts.tv_nsec = (long)((realtime_secs - rint(realtime_secs)) * 1e9);

      get_mjd_from_timespec(&ts, &(mjd->stt_imjd), &(mjd->stt_smjd), &(mjd->stt_offs));

      hputu4(st->buf, "STT_IMJD", mjd->stt_imjd);
      hputu4(st->buf, "STT_SMJD", mjd->stt_smjd);
      hputr8(st->buf, "STT_OFFS", mjd->stt_offs);
    }
    else if(state == IDLE && sttvalid != 0) {
      hputu4(st->buf, "STTVALID", 0);
    }
  }
  hashpipe_status_unlock_safe(st);
}

// Check the given block's start/stop values against the observation's OBSSTART/OBSSTOP
// values (given as PKT sequence numbers). Logic goes something like this:
//   if OBSSTART <= block_start_or_stop < OBSSTOP
//     return RECORD
//   else if block_stop <= OBSSTART
//     return ARMED
//   else
//     return IDLE
//   endif
static
enum run_states state_from_start_stop(const uint64_t obs_start_pktidx, const uint64_t obs_stop_pktidx,
                                       const uint64_t block_start_pktidx, const uint64_t block_stop_pktidx)
{
  if(block_start_pktidx >= obs_stop_pktidx){
    return IDLE;
  }
  else if(block_stop_pktidx <= obs_start_pktidx){ //block_stop_pktidx is exclusive
    return ARMED;
  }
  else{// block_start_pktidx < obs_stop_pktidx && block_stop_pktidx >= obs_start_pktidx
    return RECORD;
  }
}

int wait_out_free_transfer_in_and_set_filled(
    hpguppi_input_databuf_t* indb, int* blk_in_idx,
    hpguppi_input_databuf_t* outdb, int* blk_out_idx,
    hashpipe_status_t* st, const char* status_key
    ){
  int rv;
  // Waiting for output
  while ((rv=hpguppi_input_databuf_wait_free(outdb, *blk_out_idx)) !=
      HASHPIPE_OK)
  {
    if (rv == HASHPIPE_TIMEOUT)
    {
      hashpipe_status_lock_safe(st);
        hputs(st->buf, status_key, "outblocked");
      hashpipe_status_unlock_safe(st);
      return rv;
    }
    else
    {
      hashpipe_error(__FUNCTION__, "error waiting for output buffer, rv: %i", rv);
      return rv;
    }
  }

  memcpy(hpguppi_databuf_header(outdb, *blk_out_idx),
    hpguppi_databuf_header(indb, *blk_in_idx),
    HASHPIPE_STATUS_TOTAL_SIZE
  );
  memcpy(hpguppi_databuf_data(outdb, *blk_out_idx),
    hpguppi_databuf_data(indb, *blk_in_idx),
    BLOCK_DATA_SIZE
  );

  hpguppi_input_databuf_set_filled(outdb, *blk_out_idx);
  *blk_out_idx = (*blk_out_idx + 1) % outdb->header.n_block;

  return rv;
}

static void *run(hashpipe_thread_args_t * args)
{
    // Local aliases to shorten access to args fields
    // Our output buffer happens to be a hpguppi_input_databuf
    
  hpguppi_input_databuf_t *indb  = (hpguppi_input_databuf_t *)args->ibuf;
  hpguppi_input_databuf_t *outdb = (hpguppi_input_databuf_t *)args->obuf;

  hashpipe_status_t* st = &(args->st);
  const char* status_key = args->thread_desc->skey;
  const char* thread_name = args->thread_desc->name;

  int rv;

  int curblock_in=0;
  int curblock_out=0;
  char *datablock_header;

  struct mjd_t *mjd = malloc(sizeof(struct mjd_t));

  /* Misc counters, etc */
  uint64_t obs_npacket_total=0, obs_ndrop_total=0;
  uint32_t block_npacket=0, block_ndrop=0;

  uint64_t obs_start_pktidx = 0, obs_stop_pktidx = 0;
  uint64_t block_start_pktidx = 0, block_stop_pktidx = 0;
  
  char waiting=-1, flag_state_update=0;
  enum run_states state = IDLE;
  
  /* Heartbeat variables */
  time_t lasttime = 0;
  time_t curtime = 0;
  char timestr[32] = {0};
    
  uint32_t blocks_per_second = 0;

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
              hputu4(st->buf, "BLKSPS", blocks_per_second);
              hputs(st->buf, "DAQPULSE", timestr);
              HPUT_DAQ_STATE(st, state);
          }
          hashpipe_status_unlock_safe(st);
          blocks_per_second = 0;
      }

      // Waiting for input
      while ((rv=hpguppi_input_databuf_wait_filled(indb, curblock_in)) != 
          HASHPIPE_OK)
      {
        if (rv == HASHPIPE_TIMEOUT)
        {
          hashpipe_status_lock_safe(st);
            hputs(st->buf, status_key, "waiting");
          hashpipe_status_unlock_safe(st);
          waiting=1;
          break;
        }
        else
        {
          hashpipe_error(thread_name, "error waiting for input buffer, rv: %i", rv);
          pthread_exit(NULL);
          break;
        }
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
    hgetu8(datablock_header, "OBSSTART", &obs_start_pktidx);
    hgetu8(datablock_header, "OBSSTOP", &obs_stop_pktidx);
    hgetu8(datablock_header, "PKTSTART", &block_start_pktidx);
    hgetu8(datablock_header, "PKTSTOP", &block_stop_pktidx);

    switch(state_from_start_stop(obs_start_pktidx, obs_stop_pktidx, block_start_pktidx, block_stop_pktidx)){
      case IDLE:// If should IDLE, 
        if(state != IDLE){
          if(state == RECORD){//and recording, finalise block
            update_stt_status_keys(st, IDLE, obs_start_pktidx, mjd);
            // trigger the rawdisk thread to close the fd ("recording stopped")
            // thusly this block is not actually recorded.
            hputu8(datablock_header, "PKTIDX", block_stop_pktidx);

            wait_out_free_transfer_in_and_set_filled(indb, &curblock_in, outdb, &curblock_out, st, status_key);
          }
          flag_state_update = 1;
          state = IDLE;
        }
        break;
      case RECORD:// If should RECORD
        if (state != RECORD){
          obs_npacket_total = 0;
          obs_ndrop_total = 0;
          if(state != ARMED){// didn't arm correctly
            update_stt_status_keys(st, state, obs_start_pktidx, mjd);
            hputu4(datablock_header, "STT_IMJD", mjd->stt_imjd);
            hputu4(datablock_header, "STT_SMJD", mjd->stt_smjd);
            hputr8(datablock_header, "STT_OFFS", mjd->stt_offs);
          }
          flag_state_update = 1;
          state = RECORD;
        }
        wait_out_free_transfer_in_and_set_filled(indb, &curblock_in, outdb, &curblock_out, st, status_key);
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

    blocks_per_second ++;
    hpguppi_input_databuf_set_free(indb, curblock_in);
    curblock_in  = (curblock_in + 1) % indb->header.n_block;

    /* Will exit if thread has been cancelled */
    pthread_testcancel();
  }

  return NULL;
}

static hashpipe_thread_desc_t control_thread = {
    name: "hpguppi_atasnap_control_thread",
    skey: "CTRLSTAT",
    init: NULL,
    run:  run,
    ibuf_desc: {hpguppi_input_databuf_create},
    obuf_desc: {hpguppi_input_databuf_create}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&control_thread);
}

// vi: set ts=8 sw=4 et :
