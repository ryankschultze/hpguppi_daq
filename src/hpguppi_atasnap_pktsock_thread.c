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

#define ELAPSED_NS(start,stop) \
  (((int64_t)stop.tv_sec-start.tv_sec)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

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
enum pkt_obs_code { PKT_OBS_OK=0,
                    PKT_OBS_IDX,
                    PKT_OBS_FENG,
                    PKT_OBS_SCHAN,
                    PKT_OBS_STREAM
                  };

// These can be re-enabled once they are used
#if 0
static const uint64_t START_OK_MARGIN   =      64;
static const uint64_t START_LATE_MARGIN = (1<<20);
#endif

/* Structs/functions to more easily deal with multiple
 * active blocks being filled
 */
struct datablock_stats {
    struct hpguppi_input_databuf *dbout; // Pointer to overall shared mem databuf
    int block_idx;                    // Block index number in databuf
    int64_t block_num;                // Absolute block number
    uint64_t packet_idx;              // Index of first packet number in block
    int pktidx_per_block;            // Total number of packets to go in the block
    uint64_t pkts_per_block;
    int npacket;                      // Number of packets filled so far
    int ndrop;                     // Number of dropped packets so far
    uint64_t last_pkt;                // Last packet seq number written to block
};

// Returns pointer to datablock_stats's output data block
static char * datablock_stats_data(const struct datablock_stats *d)
{
  return hpguppi_databuf_data(d->dbout, d->block_idx);
}

// Returns pointer to datablock_stats's header
static char * datablock_stats_header(const struct datablock_stats *d)
{
  return hpguppi_databuf_header(d->dbout, d->block_idx);
}

// Reset counter(s) in datablock_stats
static void reset_datablock_stats(struct datablock_stats *d)
{
  d->npacket=0;
  d->ndrop=0;
}

// (Re-)initialize some or all fields of datablock_stats bi.
// d->db is set if dbout is non-NULL.
// d->block_idx is set if block_idx >= 0.
// d->block_num is always set and the stats are always reset.
// d->pkts_per_block is set of pkt_size > 0.
static void init_datablock_stats(struct datablock_stats *d,
    struct hpguppi_input_databuf *dbout, int block_idx, int64_t block_num,
    uint64_t pkts_per_block)
{
  if(dbout) {
    d->dbout = dbout;
  }
  if(block_idx >= 0) {
    d->block_idx = block_idx;
  }
  d->block_num = block_num;
  if(pkts_per_block > 0) {
    d->pkts_per_block = pkts_per_block;
  }
  reset_datablock_stats(d);
}

// Update block's header info and set filled status (i.e. hand-off to downstream)
static void finalize_block(struct datablock_stats *d)
{
  if(d->block_idx < 0) {
    hashpipe_error(__FUNCTION__, "datablock_stats.block_idx == %d", d->block_idx);
    pthread_exit(NULL);
  }
  char *header = datablock_stats_header(d);
  char dropstat[128];
  if(d->pkts_per_block > d->npacket) {
    d->ndrop = d->pkts_per_block - d->npacket;
  }

  sprintf(dropstat, "%d/%lu", d->ndrop, d->pkts_per_block);
  hputi4(header, "NPKT", d->npacket);
  hputi4(header, "NDROP", d->ndrop);
  hputs(header, "DROPSTAT", dropstat);
  hpguppi_input_databuf_set_filled(d->dbout, d->block_idx);
}

// Advance to next block in data buffer.  This new block will contain
// absolute block block_num.
//
// NB: The caller must wait for the new data block to be free after this
// function returns!
static void increment_block(struct datablock_stats *d, int64_t block_num)
{
  if(d->block_idx < 0) {
    hashpipe_warn(__FUNCTION__,
        "datablock_stats.block_idx == %d", d->block_idx);
  }
  if(d->dbout->header.n_block < 1) {
    hashpipe_error(__FUNCTION__,
        "datablock_stats.dbout->header.n_block == %d", d->dbout->header.n_block);
    pthread_exit(NULL);
  }

  d->block_idx = (d->block_idx + 1) % d->dbout->header.n_block;
  d->block_num = block_num;
  reset_datablock_stats(d);
}

// Wait for a datablock_stats's databuf block to be free, then copy status buffer to
// block's header and clear block's data.  Calling thread will exit on error
// (should "never" happen).  Status buffer updates made after the copy to the
// block's header will not be seen in the block's header (e.g. by downstream
// threads).  Any status buffer fields that need to be updated for correct
// downstream processing of this block must be updated BEFORE calling this
// function.  Note that some of the block's header fields will be set when the
// block is finalized (see finalize_block() for details).
static void wait_for_block_free(const struct datablock_stats * d,
    hashpipe_status_t * st, const char * status_key)
{
  int rv;
  char netstat[80] = {0};
  char netbuf_status[80];
//   int netbuf_full = hpguppi_input_databuf_total_status(d->dbout);
  //struct timespec ts_sleep = {0, 10 * 1000 * 1000}; // 10 ms
//   sprintf(netbuf_status, "%d/%d", netbuf_full, d->dbout->header.n_block);//not reporting correctly
  sprintf(netbuf_status, "%d/%d", d->block_idx, d->dbout->header.n_block);

  hashpipe_status_lock_safe(st);
  {
    hgets(st->buf, status_key, sizeof(netstat), netstat);
    hputs(st->buf, status_key, "waitfree");
    hputs(st->buf, "NETBUFST", netbuf_status);
  }
  hashpipe_status_unlock_safe(st);

  while ((rv=hpguppi_input_databuf_wait_free(d->dbout, d->block_idx))
      != HASHPIPE_OK) {
    if (rv==HASHPIPE_TIMEOUT) {
    //   netbuf_full = hpguppi_input_databuf_total_status(d->dbout);
    //   sprintf(netbuf_status, "%d/%d", netbuf_full, d->dbout->header.n_block););
      hashpipe_status_lock_safe(st);
      hputs(st->buf, status_key, "outblocked");
      hputs(st->buf, "NETBUFST", netbuf_status);
      hashpipe_status_unlock_safe(st);
    } else {
      hashpipe_error("hpguppi_atasnap_pktsock_thread",
          "error waiting for free databuf");
      pthread_exit(NULL);
    }
  }

  hashpipe_status_lock_safe(st);
  {
    hputs(st->buf, status_key, netstat);
    memcpy(datablock_stats_header(d), st->buf, HASHPIPE_STATUS_TOTAL_SIZE);
  }
  hashpipe_status_unlock_safe(st);
}

unsigned check_pkt_observability(
    const struct ata_snap_obs_info * ata_oi,
    const uint64_t pkt_idx,
    const uint64_t obs_start_pktidx,
    const uint16_t feng_id,
    const int32_t stream,
    const uint16_t pkt_schan
  )
{
  if(pkt_idx < obs_start_pktidx){
    // hashpipe_error(__FUNCTION__, "pkt_idx (%lu) < (%lu) obs_start_pktidx", pkt_idx, obs_start_pktidx);
    return PKT_OBS_IDX;
  }
  if(feng_id >= ata_oi->nants){
    hashpipe_error(__FUNCTION__, "feng_id (%u) >= (%u) ata_oi->nants", feng_id, ata_oi->nants);
    return PKT_OBS_FENG;
  }
  if(pkt_schan < ata_oi->schan){
    hashpipe_error(__FUNCTION__, "pkt_schan (%d) < (%d) ata_oi->schan", pkt_schan, ata_oi->schan);
    return PKT_OBS_SCHAN;
  }
  if(stream >= ata_oi->nstrm){
    hashpipe_error(__FUNCTION__, "stream (%d) >= (%d) ata_oi->nstrm", stream, ata_oi->nstrm);
    return PKT_OBS_STREAM;
  }
  return PKT_OBS_OK;
}

// The copy_packet_data_to_databuf() function does what it says: copies packet
// data into a data buffer. This data buffer ought to be passed on to a transposition
// thread, that produces the GUPPI RAW block layout.
//
// The data buffer block is identified by the datablock_stats structure pointed to
// by the 'd' parameter.
//
// The 'ata_oi' parameter points to the observation's obs_info data.
//
// The ata_pkt parameter points to the packet.
//
// Each packet is copied as whole, and thusly the databuf has no contiguity across
// the time samples.
// The block is filled with packets with the order:
//    PKTIDX[0 ... PIPERBLK]  (Packet-Time)  Slowest
//    FENG  [0 ... NANT]      (AntennaEnum)  
//    STREAM[0 ... NSTRM]     (Packet-Freq)  Fastest
//
// where each SNAP packet has dimensions:
//    [Slowest ------> Fastest]
//    [PKTCHAN, PKTNTIME, NPOL]
//
// The below datablock_stats_data(datablock_stats_pointer) offset is:
//  First to this pktenum (time), then to  [biggest stride]
//  first location of this FID, then to
//  first location of this stream          [smallest stride]
#define COPY_PACKET_DATA_TO_DATABUF(\
      /*const struct datablock_stats*/    datablock_stats_pointer,\
      /*const struct ata_snap_pkt*/    ata_snap_pkt_pointer,\
      /*const uint64_t*/  pkt_obs_relative_idx,\
      /*const uint16_t*/  feng_id,\
      /*const int32_t*/   stream,\
      /*const uint16_t*/  pkt_schan,\
      /*const uint32_t*/  fid_stride,\
      /*const uint32_t*/  time_stride,\
      /*const uint64_t*/  pkt_payload_size,\
      /*const uint32_t*/  pkt_ntime)\
  memcpy(datablock_stats_data(datablock_stats_pointer)+(\
        ((pkt_obs_relative_idx%datablock_stats_pointer->pktidx_per_block)/pkt_ntime) * time_stride\
            +  feng_id * fid_stride\
            +  stream * pkt_payload_size\
      ),\
      ata_snap_pkt_pointer->payload, pkt_payload_size)

/* Push all blocks down a level, losing the first one */
static void block_stack_push(struct datablock_stats *d, int nblock)
{
    int i;
    for (i=1; i<nblock; i++)
        memcpy(&d[i-1], &d[i], sizeof(struct datablock_stats));
}

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

// Check the given pktidx value against the status buffer's OBSSTART/OBSSTOP
// values (given as PKT sequence numbers). Logic goes something like this:
//   if OBSSTART <= pktidx < OBSSTOP
//     return RECORD
//   else if pktidx <= OBSSTART
//     return ARMED
//   else
//     return IDLE
//   endif
static
enum run_states state_from_start_stop(const uint64_t pktidx,
                                       const uint64_t obs_start_pktidx, const uint64_t obs_stop_pktidx)
{
  if(obs_start_pktidx <= pktidx && pktidx < obs_stop_pktidx) {
    return RECORD;
  } else if (pktidx < obs_start_pktidx && obs_start_pktidx != obs_stop_pktidx) {
    return ARMED;
  } else {// pktstop <= pktidx
    return IDLE;
  }
}

int ata_snap_obs_info_read(hashpipe_status_t *st, struct ata_snap_obs_info *obs_info)
{
  int rc = 1;//obsinfo valid

  // Get any obs info from status buffer, store values
  hashpipe_status_lock_safe(st);
  {
    // Read (no change if not present)
    hgetu4(st->buf, "FENCHAN",  &obs_info->fenchan);
    hgetu4(st->buf, "NANTS",    &obs_info->nants);
    hgetu4(st->buf, "NSTRM",    &obs_info->nstrm);
    hgetu4(st->buf, "NPOL",     &obs_info->pkt_npol);
    hgetu4(st->buf, "NBITS",    &obs_info->time_nbits);
    hgetu4(st->buf, "PKTNTIME", &obs_info->pkt_ntime);
    hgetu4(st->buf, "PKTNCHAN", &obs_info->pkt_nchan);
    hgeti4(st->buf, "SCHAN",    &obs_info->schan);
    hgetr8(st->buf, "OBSBW",    &obs_info->obs_bw);
    // If obs_info is valid
    if(ata_snap_obs_info_valid(*obs_info)) {
      ata_snap_populate_block_related_fields(BLOCK_DATA_SIZE, obs_info);
    } else {
      rc = 0;
    }
  }
  hashpipe_status_unlock_safe(st);
  return rc;
}

// This function overwrites the status values with those of the provided ata_snap_obs_info
// which primarily serves to revert any changes that were made mid-observation
int ata_snap_obs_info_write(hashpipe_status_t *st, struct ata_snap_obs_info *obs_info)
{
  int rc = 1;//obsinfo valid
  uint32_t obsnchan = 0;

  // Get any obs info from status buffer, store values
  hashpipe_status_lock_safe(st);
  {
    // If obs_info is valid
    if(ata_snap_obs_info_valid(*obs_info)) {
      obsnchan = ata_snap_obsnchan(*obs_info);
      hputs(st->buf, "OBSINFO", "VALID");
    } else {
      obsnchan = 1;
      rc = 0;
      hputs(st->buf, "OBSINFO", "INVALID");
    }

    // Write (store default/invalid values if not present)
    hputu4(st->buf, "OBSNCHAN", obsnchan);
    hputu4(st->buf, "FENCHAN",  obs_info->fenchan);
    hputu4(st->buf, "NANTS",    obs_info->nants);
    hputu4(st->buf, "NSTRM",    obs_info->nstrm);
    hputu4(st->buf, "NPOL",     obs_info->pkt_npol);
    hputu4(st->buf, "NBITS",    obs_info->time_nbits);
    hputu4(st->buf, "PKTNTIME", obs_info->pkt_ntime);
    hputu4(st->buf, "PKTNCHAN", obs_info->pkt_nchan);
    hputi4(st->buf, "SCHAN",    obs_info->schan);

    hputi4(st->buf, "BLOCSIZE", obs_info->eff_block_size);
    hputu4(st->buf, "PIPERBLK", obs_info->pktidx_per_block);
    hputu4(st->buf, "PKTSIZE",  obs_info->pkt_data_size);
  }
  hashpipe_status_unlock_safe(st);
  return rc;
}

static int init(hashpipe_thread_args_t *args)
{
  	const char * thread_name = args->thread_desc->name;

    /* Non-network essential paramaters */
    int blocsize=BLOCK_DATA_SIZE;
    int directio=1;
    int nbits=4;
    int npol=2;
    int pktnchan=256;
    int nstrm=1;
    int nant=1;
    double obsbw=187.5;
    int obsnchan=64;
    int obsschan=0;
    int overlap=0;
    double tbin=0.0;
    char obs_mode[80] = {0};
    
    struct hpguppi_pktsock_params *p_psp = (struct hpguppi_pktsock_params *)
        malloc(sizeof(struct hpguppi_pktsock_params));

    if(!p_psp) {
        perror(__FUNCTION__);
        return -1;
    }

    strcpy(obs_mode, "RAW");

    hashpipe_status_t *st = &args->st;

    hashpipe_status_lock_safe(st);
    {
      // Get network parameters (BINDHOST, BINDPORT, PKTFMT)
      hpguppi_read_pktsock_params(st->buf, p_psp);
      // Get info from status buffer if present (no change if not present)
      hgeti4(st->buf, "BLOCSIZE", &blocsize);
      hgeti4(st->buf, "DIRECTIO", &directio);
      hgeti4(st->buf, "NBITS", &nbits);
      hgeti4(st->buf, "NPOL", &npol);
      hgeti4(st->buf, "PKTNCHAN", &pktnchan);
      hgeti4(st->buf, "NSTRM", &nstrm);
      hgeti4(st->buf, "NANTS", &nant);
      hgetr8(st->buf, "OBSBW", &obsbw);
      hgeti4(st->buf, "OBSSCHAN", &obsschan);
      hgeti4(st->buf, "OVERLAP", &overlap);
      hgets(st->buf, "OBS_MODE", sizeof(obs_mode), obs_mode);

      // Clean any inappropriate values 
      pktnchan = (pktnchan == 0 ? 256 : pktnchan);
      nstrm = (nstrm == 0 ? 1 : nstrm);
      nant = (nant == 0 ? 1 : nant);
      obsnchan = pktnchan*nstrm*nant;

      // Calculate TBIN
      tbin = obsnchan / obsbw / 1e6;
      // Store bind host/port info etc in status buffer
      hputs(st->buf, "BINDHOST", p_psp->ifname);
      hputi4(st->buf, "BINDPORT", p_psp->port);
      hputs(st->buf, "PKTFMT", p_psp->packet_format);
      hputi4(st->buf, "BLOCSIZE", blocsize);
      hputi4(st->buf, "DIRECTIO", directio);
      hputi4(st->buf, "NBITS", nbits);
      hputi4(st->buf, "NPOL", npol);
      hputi4(st->buf, "PKTNCHAN", pktnchan);
      hputi4(st->buf, "NSTRM", nstrm);
      hputi4(st->buf, "NANTS", nant);
      hputr8(st->buf, "OBSBW", obsbw);
      hputi4(st->buf, "OBSNCHAN", obsnchan);
      hputi4(st->buf, "OBSSCHAN", obsschan);
      hputi4(st->buf, "OVERLAP", overlap);
      hputr8(st->buf, "TBIN", tbin);
      hputs(st->buf, "OBS_MODE", obs_mode);
      // Data are in time-major order (i.e. time dimension changes faster than
      // channel dimension), so specify that data are NOT in channel major order.
      hputi4(st->buf, "CHANMAJ", 0);
    }
    hashpipe_status_unlock_safe(st);

    // Set up pktsock.  Make frame_size be a divisor of block size so that
    // frames will be contiguous in mapped mempory.  block_size must also be a
    // multiple of page_size.  Easiest way is to oversize the frames to be
    // 16384 bytes, which is bigger than we need, but keeps things easy.
    p_psp->ps.frame_size = PKTSOCK_BYTES_PER_FRAME;
    // total number of frames
    p_psp->ps.nframes = PKTSOCK_NFRAMES;
    // number of blocks
    p_psp->ps.nblocks = PKTSOCK_NBLOCKS;

    if (hashpipe_pktsock_open(&p_psp->ps, p_psp->ifname, PACKET_RX_RING) != HASHPIPE_OK) {
        hashpipe_error(thread_name, "Error opening pktsock.");
        pthread_exit(NULL);
    }

    // Store hpguppi_pktsock_params pointer in args
    args->user_data = p_psp;

    // Success!
    return 0;
}

static void *run(hashpipe_thread_args_t * args)
{
    // Local aliases to shorten access to args fields
    // Our output buffer happens to be a hpguppi_input_databuf
    hpguppi_input_databuf_t *db = (hpguppi_input_databuf_t *)args->obuf;
    hashpipe_status_t *st = &args->st;
  	const char * thread_name = args->thread_desc->name;
    const char * status_key = args->thread_desc->skey;
    struct hpguppi_pktsock_params *p_ps_params =
        (struct hpguppi_pktsock_params *)args->user_data;

    /* Read in general parameters */
    struct hpguppi_params gp;
    struct psrfits pf;
    pf.sub.dat_freqs = NULL;
    pf.sub.dat_weights = NULL;
    pf.sub.dat_offsets = NULL;
    pf.sub.dat_scales = NULL;
    char status_buf[HASHPIPE_STATUS_TOTAL_SIZE];
    hashpipe_status_lock_safe(st);
      memcpy(status_buf, st->buf, HASHPIPE_STATUS_TOTAL_SIZE);
    hashpipe_status_unlock_safe(st);
    hpguppi_read_obs_params(status_buf, &gp, &pf);
    // Structure to hold observation info, init all fields to invalid values
    pthread_cleanup_push((void *)hpguppi_free_psrfits, &pf);

    struct ata_snap_obs_info obs_info;
    ata_snap_obs_info_init(&obs_info);

    struct mjd_t *mjd = malloc(sizeof(struct mjd_t));

    int block_size = BLOCK_DATA_SIZE;

    if (hgeti4(status_buf, "BLOCSIZE", &block_size)==0) {
            block_size = BLOCK_DATA_SIZE;
            hputi4(status_buf, "BLOCSIZE", block_size);
    } else {
        if (block_size > BLOCK_DATA_SIZE) {
            hashpipe_error(thread_name, "BLOCSIZE > databuf block_size (%d > %d)", block_size, BLOCK_DATA_SIZE);
            block_size = BLOCK_DATA_SIZE;
            hputi4(status_buf, "BLOCSIZE", block_size);
        }
    }
    
    /* Figure out size of data in each packet, number of packets
     * per block, etc.  Changing packet size during an obs is not
     * recommended, so observational flexibility is only possible 
     * when hashpipe is idling.
     */
    ata_snap_obs_info_read(st, &obs_info);
    ata_snap_obs_info_write(st, &obs_info);
    
    fprintf(stderr, "Packets per block %d, Packet timestamps per block %d\n", obs_info.pkt_per_block, obs_info.pktidx_per_block);
    unsigned long pkt_blk_num, last_pkt_blk_num = ~0;


    // The incoming packets are taken from blocks of the input databuf and then
    // converted to GUPPI RAW format in blocks of the output databuf to pass to
    // the downstream thread.  We currently support two active output blocks (aka
    // "working blocks").  Working blocks are associated with absolute output
    // block numbers, which are simply PKTIDX values divided by the number of
    // packets per block (discarding any remainder).  Let the block numbers for
    // the first working block (wblk[0]) be W.  The block number for the second
    // working block (wblk[1]) will be W+1.  Incoming packets corresponding to
    // block W or W+1 are placed in the corresponding data buffer block.
    // Incoming packets for block W+2 cause block W to be "finalized" and handed
    // off to the downstream thread, working block 1 moves to working block 0 and
    // working block 1 is incremented to be W+2.  Things get "interesting" when a
    // packet is recevied for block < W or block > W+2.  Packets for block W-1
    // are ignored.  Packets with PKTIDX P corresponding block < W-1 or block >
    // W+2 cause the current working blocks' block numbers to be reset such that
    // W will refer to the block containing P and W+1 will refer to the block
    // after that.
    //
    // wblk is a two element array of datablock_stats structures (i.e. the working
    // blocks)
    /* List of databuf blocks currently in use */
    int wblk_idx;
    const int n_wblock = 2;
    struct datablock_stats wblk[n_wblock];
    // Initialize working blocks
    for(wblk_idx=0; wblk_idx<n_wblock; wblk_idx++) {
        // important to initialise the blocks with unique block idx
        // and to set the block_num < -1 so the first observation packet causes a 
        // discontinuity and so a contemporary re-init of the working blocks
        init_datablock_stats(wblk+wblk_idx, db, wblk_idx, -10, obs_info.pkt_per_block);
        wait_for_block_free(wblk+wblk_idx, st, status_key);
    }
    char *datablock_header;

    /* Misc counters, etc */
    uint64_t npacket_total=0, ndrop_total=0, nbogus_total=0, lastnpacket_total=0;
    uint64_t obs_npacket_total=0, obs_ndrop_total=0, obs_block_discontinuities=0;

    uint64_t pkt_seq_num, first_pkt_seq_num=0;
    uint64_t obs_start_pktidx = 0, obs_stop_pktidx = 0, pkt_obs_relative_idx=0;
    uint64_t pkt_payload_size;
    uint32_t fid_stride, time_stride;
    int32_t stream;
    uint16_t feng_id, pkt_schan;

    char waiting=-1;
    enum run_states state = IDLE;
    uint32_t subsequent_state_idle_count = 0;
    char flag_state_update = 0;
    char  PKT_OBS_IDX_flagged, PKT_OBS_FENG_flagged, 
          PKT_OBS_SCHAN_flagged, PKT_OBS_STREAM_flagged, pkt_obs_code;
    char flag_obs_end = 0, flag_arm_for_obs = 0;

    /* Time parameters */
    // int stt_imjd=0, stt_smjd=0;
    // double stt_offs=0.0;
    struct timespec ts_start_block = {0}, ts_stop_block = {0};
    struct timespec ts_checked_obs_startstop = {0}, ts_now = {0};
    // Heartbeat variables
    time_t lasttime = 0;
    time_t curtime = 0;
    char timestr[32] = {0};
    
    float blocks_per_second = 0.0;

    // Drop all packets to date
    unsigned char *p_frame;
    while((p_frame=hashpipe_pktsock_recv_frame_nonblock(&p_ps_params->ps))) {
        hashpipe_pktsock_release_frame(p_frame);
    }
    struct ata_snap_pkt *ata_snap_pkt;

    fprintf(stderr, "Receiving at interface %s, port %d\n",
    p_ps_params->ifname, p_ps_params->port);//p_ps_params->packet_size);

    /* Main loop */
    while (run_threads()) {

        /* Wait for data */
        do {
            // Heartbeat update?
            time(&curtime);//time stores seconds since epoch
            if(flag_state_update || curtime > lasttime) {// once per second
                if (state == IDLE){
                  ata_snap_obs_info_read(st, &obs_info);
                }
                ata_snap_obs_info_write(st, &obs_info);
                flag_state_update = 0;
                lasttime = curtime;

                ctime_r(&curtime, timestr);
                timestr[strlen(timestr)-1] = '\0'; // Chop off trailing newline
                hashpipe_status_lock_safe(st);
                {
                    hputi8(st->buf, "NPKTS", npacket_total);
                    hputi8(st->buf, "OBSDSCNT", obs_block_discontinuities);
                    hputi8(st->buf, "OBSNPKTS", obs_npacket_total);
                    hputi8(st->buf, "OBSNDROP", obs_ndrop_total);
                    hputi8(st->buf, "NDROP", ndrop_total);
                    hputr4(st->buf, "BLKSPS", blocks_per_second);
                    hputr4(st->buf, "PHYSPKPS", npacket_total-lastnpacket_total);
                    hputr4(st->buf, "PHYSGBPS", ((npacket_total-lastnpacket_total)*obs_info.pkt_data_size)/1e9);
                    hputs(st->buf, "DAQPULSE", timestr);
                    HPUT_DAQ_STATE(st, state);

                    // Overwrite any changes to reflect the static nature of these keys
                    hputs(st->buf, "BINDHOST", p_ps_params->ifname);
                    hputi4(st->buf, "BINDPORT", p_ps_params->port);
                }
                hashpipe_status_unlock_safe(st);
                lastnpacket_total = npacket_total;
            }

            p_frame = hashpipe_pktsock_recv_udp_frame(
                &p_ps_params->ps, p_ps_params->port, 500); // 0.5 second timeout

            /* Set "waiting" flag */
            if (!p_frame && run_threads() && waiting!=1) {
                hashpipe_status_lock_safe(st);
                hputs(st->buf, status_key, "waiting");
                hashpipe_status_unlock_safe(st);
                waiting=1;
            }

        } while (!p_frame && run_threads());

        if(!run_threads()) {
            // We're outta here!
            hashpipe_pktsock_release_frame(p_frame);
            break;
        }

        /* Check packet size */
        if(obs_info.pkt_data_size != PKT_UDP_SIZE(p_frame) - 8) {
            /* Unexpected packet size, ignore */
            hashpipe_status_lock_safe(st);
              hputi4(st->buf, "NBOGUS", ++nbogus_total);
              hputi4(st->buf, "BOGUSIZE", PKT_UDP_SIZE(p_frame)-8);
            hashpipe_status_unlock_safe(st);
            // Release frame!
            hashpipe_pktsock_release_frame(p_frame);
            continue;
        }
        
        /* Update status if needed */
        if (waiting!=0) {
            hashpipe_status_lock_safe(st);
            hputs(st->buf, status_key, "processing");
            hashpipe_status_unlock_safe(st);
            waiting=0;
        }

        ata_snap_pkt = (struct ata_snap_pkt*) p_frame;
        // Get packet's sequence number
        pkt_seq_num =  ATA_SNAP_PKT_NUMBER(ata_snap_pkt);

        // Catch any changes in OBSSTART/STOP
        clock_gettime(CLOCK_MONOTONIC, &ts_now);
        if(ELAPSED_NS(ts_checked_obs_startstop, ts_now) > 50*1000*1000){
          memcpy(&ts_checked_obs_startstop, &ts_now, sizeof(struct timespec));
          hashpipe_status_lock_safe(st);
          {
            hgetu8(st->buf, "OBSSTART", &obs_start_pktidx);
            hgetu8(st->buf, "OBSSTOP", &obs_stop_pktidx);
          }
          hashpipe_status_unlock_safe(st);
        }
          
        switch(state_from_start_stop(pkt_seq_num, obs_start_pktidx, obs_stop_pktidx)){
          case IDLE:// If should IDLE, 
            flag_state_update = (state != IDLE ? 1 : 0);// flag state update
            subsequent_state_idle_count += (state == RECORD ? 1 : 0);
            if(state == RECORD && subsequent_state_idle_count > 100){//and recording, finalise block
              flag_obs_end = 1;
              // first_pkt_seq_num = 0; // reset after finalisation of block
              state = IDLE;
              update_stt_status_keys(st, state, pkt_seq_num, mjd);
              // update PKTIDX triggering rawdisk to close fd when this last block gets finalised.
              hputi8(datablock_stats_header(wblk), "PKTIDX", pkt_seq_num);
            }
            else if(state == ARMED){
              state = IDLE;
            }
            break;
          case RECORD:// If should RECORD
            subsequent_state_idle_count = 0;
            if (state != RECORD && ata_snap_obs_info_valid(obs_info)){// Only enter recording mode if obs_params are valid
              // and not recording, flag obs_start
              flag_state_update = 1;
              flag_arm_for_obs = (state == IDLE ? 1 : 0);
              state = RECORD;
            }
            break;
          case ARMED:// If should ARM,
            subsequent_state_idle_count = 0;
            if(state == IDLE){
              flag_state_update = 1;// flag state update
              state = ARMED;
              flag_arm_for_obs = 1;
            }
          default:
            break;
        }
        
        if(flag_arm_for_obs == 1){
          flag_arm_for_obs = 0;

          PKT_OBS_IDX_flagged = 0;
          PKT_OBS_FENG_flagged = 0;
          PKT_OBS_SCHAN_flagged = 0;
          PKT_OBS_STREAM_flagged = 0;
          
          first_pkt_seq_num = obs_start_pktidx;
          ata_snap_obs_info_read(st, &obs_info);
          obs_npacket_total = 0;
          obs_ndrop_total = 0;
          obs_block_discontinuities = 0;
          update_stt_status_keys(st, state, first_pkt_seq_num, mjd);

          pkt_payload_size = ata_snap_pkt_payload_bytes(obs_info);
          fid_stride = obs_info.nstrm*pkt_payload_size;
          time_stride = obs_info.nants*fid_stride;

          for(wblk_idx=0; wblk_idx<n_wblock; wblk_idx++) {
            wblk[wblk_idx].pkts_per_block = obs_info.pkt_per_block;
            wblk[wblk_idx].pktidx_per_block = obs_info.pktidx_per_block;
            init_datablock_stats(wblk+wblk_idx, NULL, -1,
                0+wblk_idx,
                obs_info.pkt_per_block);

            // also update the working blocks' headers
            datablock_header = datablock_stats_header(wblk+wblk_idx);
            hashpipe_status_lock_safe(st);
            {
              memcpy(datablock_header, st->buf, HASHPIPE_STATUS_TOTAL_SIZE);
            }
            hashpipe_status_unlock_safe(st);
            hputu8(datablock_header, "PKTIDX", first_pkt_seq_num + wblk_idx * obs_info.pktidx_per_block);
            hputu8(datablock_header, "PKTSTART", first_pkt_seq_num + wblk_idx * obs_info.pktidx_per_block);
            hputu8(datablock_header, "PKTSTOP", first_pkt_seq_num + (wblk_idx + 1) * obs_info.pktidx_per_block);
          }
          hashpipe_info(thread_name, "Armed wblk for observation: first_pkt_seq_num = %lu", first_pkt_seq_num);
        }

        pkt_obs_relative_idx = pkt_seq_num - first_pkt_seq_num;
        pkt_blk_num = pkt_obs_relative_idx / obs_info.pktidx_per_block;
        // fprintf(stderr, "seq: %012ld\tblk: %06ld\n", pkt_seq_num, pkt_blk_num);

        // Update PKTIDX in status buffer if pkt_seq_num % obs_info.pktidx_per_block == 0
        if(pkt_blk_num != last_pkt_blk_num){

          last_pkt_blk_num = pkt_blk_num;

          hashpipe_status_lock_safe(st);
          hputu8(st->buf, "BLKIDX", pkt_blk_num);
          hputu8(st->buf, "PKTIDX", pkt_seq_num);
          hputu8(st->buf, "PKTSTART", pkt_seq_num);
          hputu8(st->buf, "PKTSTOP", pkt_seq_num + obs_info.pktidx_per_block);
          hashpipe_status_unlock_safe(st);
        }

        // Tally npacket_totals
        npacket_total += 1;
        // count ndrop only when a block is finalised, lest it is filled out of order.
        
        if (state != RECORD && flag_obs_end != 1){
          hashpipe_pktsock_release_frame(p_frame);
          continue;
        }

        feng_id = ATA_SNAP_PKT_FENG_ID(ata_snap_pkt);
        pkt_schan = ATA_SNAP_PKT_CHAN(ata_snap_pkt);
        stream = (pkt_schan - obs_info.schan) / obs_info.pkt_nchan;
        pkt_obs_code = check_pkt_observability(&obs_info,// ata_snap_pkt, 
                    pkt_seq_num, first_pkt_seq_num, feng_id, stream, pkt_schan);

        // Manage blocks based on pkt_blk_num
        // fprintf(stderr, "%010ld\r", pkt_blk_num);
        if(flag_obs_end || pkt_blk_num == wblk[n_wblock-1].block_num + 1) {
            // Time to advance the blocks!!!
            // fprintf(stderr, "\nFinalising Block: %ld", wblk[0].block_num);
            clock_gettime(CLOCK_MONOTONIC, &ts_stop_block);
            // Finalize first working block
            datablock_header = datablock_stats_header(&wblk[0]);
            hputu8(datablock_header, "PKTIDX", first_pkt_seq_num + wblk[0].block_num * obs_info.pktidx_per_block);
            hputu8(datablock_header, "PKTSTART", first_pkt_seq_num + wblk[0].block_num * obs_info.pktidx_per_block);
            hputu8(datablock_header, "PKTSTOP", first_pkt_seq_num + (wblk[0].block_num + 1) * obs_info.pktidx_per_block);
            finalize_block(wblk);
            if(!flag_obs_end){
              // Update ndrop counter
              obs_ndrop_total += wblk->ndrop;
              obs_npacket_total += wblk->npacket;
              ndrop_total += wblk->ndrop;
            }
            // Shift working blocks
            block_stack_push(wblk, n_wblock);
            // Increment last working block
            increment_block(&wblk[n_wblock-1], pkt_blk_num);
            // Wait for new databuf data block to be free
            wait_for_block_free(&wblk[n_wblock-1], st, status_key);
            
            blocks_per_second = 1000.0*1000.0*1000.0/ELAPSED_NS(ts_start_block,ts_stop_block);
            clock_gettime(CLOCK_MONOTONIC, &ts_start_block);
        }
        // Check for PKTIDX discontinuity
        else if(pkt_obs_code == PKT_OBS_OK && 
                  (pkt_blk_num + 1 < wblk[0].block_num
                || pkt_blk_num > wblk[n_wblock-1].block_num + 1)
                ) {
            // Should only happen when transitioning into ARMED, so warn about it
            hashpipe_warn(thread_name,
                "working blocks reinit due to packet discontinuity\n\t\t(PKTIDX %lu) [%ld, %ld  <> %lu]",
                pkt_seq_num, wblk[0].block_num - 1, wblk[n_wblock-1].block_num + 1, pkt_blk_num);
            obs_block_discontinuities++;
            // Re-init working blocks for block number of current packet's block,
            // and clear their data buffers
            for(wblk_idx=0; wblk_idx<n_wblock; wblk_idx++) {
              init_datablock_stats(wblk+wblk_idx, NULL, -1,
                  pkt_blk_num+wblk_idx,
                  obs_info.pkt_per_block);

              // also update the working blocks' headers
              datablock_header = datablock_stats_header(wblk+wblk_idx);
              memcpy(datablock_header, st->buf, HASHPIPE_STATUS_TOTAL_SIZE);
              hputu8(datablock_header, "PKTIDX", pkt_seq_num + wblk_idx * obs_info.pktidx_per_block);
              hputu8(datablock_header, "PKTSTART", pkt_seq_num + wblk_idx * obs_info.pktidx_per_block);
              hputu8(datablock_header, "PKTSTOP", pkt_seq_num + (wblk_idx + 1) * obs_info.pktidx_per_block);
            }
            // Immediately update STT keys to ensure that the rawdisk thread uses a new stem
            update_stt_status_keys(st, state, obs_start_pktidx, mjd);
        }

        // Check observation state
        if(state != RECORD){// Either ARMED or transitioning RECORD->IDLE
          if(flag_obs_end){
            flag_obs_end = 0;
            first_pkt_seq_num = 0;
          }
          hashpipe_pktsock_release_frame(p_frame);
          continue;
        }

        // Only copy packet data and count packet if its wblk_idx is valid
        switch(pkt_obs_code){
          case PKT_OBS_OK:
            // Once we get here, compute the index of the working block corresponding
            // to this packet.  The computed index may not correspond to a valid
            // working block!
            wblk_idx = pkt_blk_num - wblk[0].block_num;

            if(0 <= wblk_idx && wblk_idx < n_wblock) {
              // Copy packet data to data buffer of working block
              COPY_PACKET_DATA_TO_DATABUF(((struct datablock_stats*) wblk+wblk_idx),
                  ata_snap_pkt, pkt_obs_relative_idx,
                  feng_id, stream, pkt_schan,
                  fid_stride, time_stride, pkt_payload_size, 16);//obs_info.pkt_ntime);
              // Count packet for block and for processing stats
              wblk[wblk_idx].npacket++;
            }
            else{
              hashpipe_error(thread_name, "Packet ignored: determined wblk_idx = %d", wblk_idx);
            }
            break;
          case PKT_OBS_IDX:
            if(!PKT_OBS_IDX_flagged){
              PKT_OBS_IDX_flagged = 1;
              hashpipe_error(thread_name, "Packet ignored: PKT_OBS_IDX %d", pkt_seq_num);
            }
            break;
          case PKT_OBS_FENG:
            if(!PKT_OBS_FENG_flagged){
              PKT_OBS_FENG_flagged = 1;
              hashpipe_error(thread_name, "Packet ignored: PKT_OBS_FENG");
            }
            break;
          case PKT_OBS_SCHAN:
            if(!PKT_OBS_SCHAN_flagged){
              PKT_OBS_SCHAN_flagged = 1;
              hashpipe_error(thread_name, "Packet ignored: PKT_OBS_SCHAN");
            }
            hashpipe_status_lock_safe(st);
            hputs(st->buf, "OBSINFO", "INVALID SCHAN");
            hashpipe_status_unlock_safe(st);
            break;
          case PKT_OBS_STREAM:
            if(!PKT_OBS_STREAM_flagged){
              PKT_OBS_STREAM_flagged = 1;
              hashpipe_error(thread_name, "Packet ignored: PKT_OBS_STREAM");
            }
            hashpipe_status_lock_safe(st);
            hputs(st->buf, "OBSINFO", "INVALID NSTRM");
            hashpipe_status_unlock_safe(st);
            break;
          default:
            break;
        }

        // Release frame back to ring buffer
        hashpipe_pktsock_release_frame(p_frame);

        /* Will exit if thread has been cancelled */
        pthread_testcancel();
    }

    pthread_exit(NULL);

    /* Have to close all push's */
    pthread_cleanup_pop(0); /* Closes hpguppi_free_psrfits */

    return NULL;
}

static hashpipe_thread_desc_t net_thread = {
    name: "hpguppi_atasnap_pktsock_thread",
    skey: "NETSTAT",
    init: init,
    run:  run,
    ibuf_desc: {NULL},
    obuf_desc: {hpguppi_input_databuf_create}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&net_thread);
}




#if 0 // this transposition has been outsourced to pkt_to_FTP_transpose
// The copy_packet_data_timecontiguously_to_databuf() function does what it says: copies packet
// data into a data buffer.
//
// The data buffer block is identified by the datablock_stats structure pointed to
// by the bi parameter.
//
// The p_oi parameter points to the observation's obs_info data.
//
// The p_fei parameter points to an ata_snap_feng_info structure containing the
// packet's metadata.
//
// The p_payload parameter points to the payload of the packet.
//
// This is for packets in [channel (slowest), time, pol (fastest)] order.  In
// other words:
//     t=0           t=1               t=pkt_ntime-1
//     C0T0P0 C0T0P1 C0T1P0 C0T1P1 ... C0TtP0 C0TtP1 <- c = 0
//     C1T0P0 C1T0P1 C1T1P0 C1T1P1 ... C1TtP0 C1TtP1 <- c = 1
//     ...
//     CcT0P0 CcT0P1 CcT1P0 CcT1P1 ... CcTtP0 CcTtP1 <- c = pkt_nchan
//
// GUPPI RAW block is ordered as:
//
//     t=0               t=1                   t=pktidx_per_block
//     F0C0T0P0 F0C0T0P1 F0C0T1P0 F0C0T1P1 ... F0C0TtP0 F0C0TtP1
//     F0C1T1P0 F0C1T1P1 F0C1T1P0 F0C1T1P1 ... F0C1TtP0 F0C1TtP1
//     ...
//     F0CcT0P0 F0CcT0P1 F0CcT1P0 F0CcT1P1 ... F0CcTtP0 F0CcTtP1
//     F1C0T0P0 F1C0T0P1 F1C0T1P0 F1C0T1P1 ... F1C0TtP0 F1C0TtP1
//     F1C1T1P0 F1C1T1P1 F1C1T1P0 F1C1T1P1 ... F1C1TtP0 F1C1TtP1
//     ...
//     ...
//     FfCcT0P0 FfCcT0P1 FfC0T1P0 FfCcT1P1 ... FfCcTtP0 FfCcTtP1
//
// where F is FID (f=NANTS-1), T is time (t=pktidx_per_block-1), C is channel
// (c=NSTRMS*PKT_NCHAN-1), and P is polarization.  Streams are not shown
// separately, which is why they are bundled in with channel number.  Each
// packet's channels need to be split out over the GUPPI RAW block, which can be
// shown somewhat pictorially like this (with time faster changing in the horizontal
// direction, then channel changing in the vertical direction) for a single
// PKTIDX value (i.e. this is a slice in time of the GUPPI RAW block):
//
//     [FID=0, STREAM=0, CHAN=0:PKT_NCHAN-1, TIME=0:PKT_NTIME-1] ...
//     [FID=0, STREAM=1, CHAN=0:PKT_NCHAN-1, TIME=0:PKT_NTIME-1] ...
//      ...
//     [FID=0, STREAM=s, CHAN=0:PKT_NCHAN-1, TIME=0:PKT_NTIME-1] ...
//     [FID=1, STREAM=0, CHAN=0:PKT_NCHAN-1, TIME=0:PKT_NTIME-1] ...
//     [FID=1, STREAM=1, CHAN=0:PKT_NCHAN-1, TIME=0:PKT_NTIME-1] ...
//      ...
//      ...
//     [FID=f, STREAM=s, CHAN=0:PKT_NCHAN-1, TIME=0:PKT_NTIME-1] ...
//
static char copy_packet_printed = 0;
static void copy_packet_data_timecontiguously_to_databuf(const struct datablock_stats *d,
    const struct ata_snap_obs_info * ata_oi,
    struct ata_snap_pkt* ata_pkt, const uint64_t obs_start_pktidx)
{
    // the pointer's data width means the offset is in terms of bytes
    char *dst_base = datablock_stats_data(d);
    unsigned char *pkt_payload = ata_pkt->payload;
    const uint16_t pkt_schan = ATA_SNAP_PKT_CHAN(ata_pkt);
    const uint16_t feng_id = ATA_SNAP_PKT_FENG_ID(ata_pkt);
    // offset pkt_idx so it the first of the observation is at the start of the block
    const uint64_t pkt_idx =  ATA_SNAP_PKT_NUMBER(ata_pkt) - obs_start_pktidx;

    // The packet's individual channels have a data size partial to the packet's payload size.
    // The pkt_channel_size includes all time samples and polarisations for the channel.
    const uint32_t pkt_channel_size = ata_snap_pkt_payload_bytes(*ata_oi)/ata_oi->pkt_nchan;

    // blk_channel_stride is the length of a channel's data within the RAW data block,
    // which equals the number of packet indices per block, as the packet index increments
    // in steps of pkt_ntime.
    const uint32_t blk_channel_stride = d->pktidx_per_block*ata_oi->pkt_npol;

    // stream_stride is the size of pkt_nchan channels
    const uint32_t stream_stride = blk_channel_stride*ata_oi->pkt_nchan;

    // fid_stride is the size of all streams of a single F engine:
    const uint32_t fid_stride = stream_stride * ata_oi->nstrm;

    // Stream is the channel grouping for an antenna
    const uint32_t stream = (pkt_schan - ata_oi->schan) / ata_oi->pkt_nchan;

    // Advance dst_base to... 
    const uint64_t pktchan0_offset = feng_id * fid_stride // first location of this FID, then
            +  stream * stream_stride                     // first location of this stream
            +  ((pkt_idx%d->pktidx_per_block)*ata_oi->pkt_npol);
    dst_base += pktchan0_offset;

    if(!copy_packet_printed || pktchan0_offset < 0 || pktchan0_offset > 128*1024*1024){
      printf("pkt_channel_size = %d\n", pkt_channel_size);
      printf("pkt_nchan        = %d\n", ata_oi->pkt_nchan);
      printf("pkt_schan        = %d\n", pkt_schan);
      printf("schan            = %d\n", ata_oi->schan);
      printf("stream           = %d\n", stream);
      printf("feng_id          = %d\n", feng_id);
      printf("stream_stride    = %u\n", stream_stride);
      printf("fid_stride       = %u\n", fid_stride);
      printf("pktchan0_offset  = %lu\n", pktchan0_offset);
      if (!copy_packet_printed && pkt_idx != 0){
        copy_packet_printed = 1;
      }
      else if(pkt_idx != 0){
        exit(1);
      }
    }

    for(int pkt_chan = 0; pkt_chan < ata_oi->pkt_nchan; pkt_chan++){
      memcpy(dst_base, pkt_payload, pkt_channel_size);
      dst_base += blk_channel_stride;
      pkt_payload += pkt_channel_size;
    }

    return;
}
#endif

// vi: set ts=8 sw=4 et :
