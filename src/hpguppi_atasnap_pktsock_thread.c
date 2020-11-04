/* hpguppi_net_thread.c
 *
 * Routine to read packets from network and put them
 * into shared memory blocks.
 */

#define _GNU_SOURCE 1
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

#define PKTS_PER_SECOND_BUF_LENGTH 5

#define PKTSOCK_BYTES_PER_FRAME (16384)
#define PKTSOCK_FRAMES_PER_BLOCK (8)
#define PKTSOCK_NBLOCKS (800)
#define PKTSOCK_NFRAMES (PKTSOCK_FRAMES_PER_BLOCK * PKTSOCK_NBLOCKS)

#define ELAPSED_NS(start,stop) \
  (((int64_t)stop.tv_sec-start.tv_sec)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

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

// These can be re-enabled once they are used
#if 0
static const uint64_t START_OK_MARGIN   =      64;
static const uint64_t START_LATE_MARGIN = (1<<20);
#endif

/* It's easier to just make these global ... */
static uint64_t npacket_total=0, ndrop_total=0, nbogus_total=0;

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
  hputi8(header, "PKTIDX", d->block_num * d->pktidx_per_block);
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
      hashpipe_error("hpguppi_atasnap_voltage_thread",
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

// The copy_packet_data_to_databuf() function does what it says: copies packet
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
// This is for packets in [time (slowest), channel, pol (fastest)] order.  In
// other words:
//
//     T0C0P0 T0C0P1 T0C1P0 T0C1P1 ... T0CcP0 T0CcP1 <- t=0
//     T1C0P0 T1C0P1 T1C1P0 T1C1P1 ... T1CcP0 T1CcP1 <- t=1
//     ...
//     TtC0P0 TtC0P1 TtC1P0 TtC1P1 ... TtCcP0 TtCcP1 <- t=pkt_ntime-1
//
// GUPPI RAW block is ordered as:
//
//     t=0               t=1                   t=NTIME
//     F0T0C0P0 F0T0C0P1 F0T1C0P0 F0T1C0P1 ... F0TtC0P0 F0TtC0P1
//     F0T0C1P0 F0T0C1P1 F0T1C1P0 F0T1C1P1 ... F0TtC1P0 F0TtC1P1
//     ...
//     F0T0CcP0 F0T0CcP1 F0T1CcP0 F0T1CcP1 ... F0TtCcP0 F0TtCcP1
//     F1T0C0P0 F1T0C0P1 F1T1C0P0 F1T1C0P1 ... F1TtC0P0 F1TtC0P1
//     F1T0C1P0 F1T0C1P1 F1T1C1P0 F1T1C1P1 ... F1TtC1P0 F1TtC1P1
//     ...
//     ...
//     FfT0CcP0 FfT0CcP1 FfT1CcP0 FfT1CcP1 ... FfTtCcP0 FfTtCcP1
//
// where F is FID (f=NANTS-1), T is time (t=PKT_NTIME-1), C is channel
// (c=NSTRMS*PKT_NCHAN-1), and P is polarization.  Streams are not shown
// separately, which is why they are bundled in with channel number.  Each
// packet fills a 2D rectangle in the GUPPI RAW block, which can be shown
// somewhat pictorially like this (with time faster changing in the horizontal
// direction, then channel changing in the vertical direction) for a single
// PKTIDX value (i.e. this is a slice in time of the GUPPI RAW block):
//
//     [FID=0, STREAM=0, TIME=0:PKT_NTIME-1, CHAN=0:PKT_NCHAN-1] ...
//     [FID=0, STREAM=1, TIME=0:PKT_NTIME-1, CHAN=0:PKT_NCHAN-1] ...
//      ...
//     [FID=0, STREAM=s, TIME=0:PKT_NTIME-1, CHAN=0:PKT_NCHAN-1] ...
//     [FID=1, STREAM=0, TIME=0:PKT_NTIME-1, CHAN=0:PKT_NCHAN-1] ...
//     [FID=1, STREAM=1, TIME=0:PKT_NTIME-1, CHAN=0:PKT_NCHAN-1] ...
//      ...
//      ...
//     [FID=f, STREAM=s, TIME=0:PKT_NTIME-1, CHAN=0:PKT_NCHAN-1] ...
//
// int copy_packet_data_to_databuf_printed = 1;
static void copy_packet_data_to_databuf(const struct datablock_stats *d,
    const struct ata_snap_obs_info * ata_oi,
    struct ata_snap_pkt* ata_pkt, int pkt_payload_size)
{
    uint16_t * dst_base = (uint16_t *)datablock_stats_data(d);
    unsigned long pkt_idx = __bswap_64(ata_pkt->timestamp);
    unsigned short feng_id = __bswap_16(ata_pkt->feng_id);

    // stream_stride is the size of a single stream for a single F engine for all
    // NTIME samples of the block and all channels in a stream (i.e. in a packet):
    const int stream_stride = pkt_payload_size * d->pktidx_per_block;

    // fid_stride is the size of all streams of a single F engine:
    const int fid_stride = stream_stride * ata_oi->nstrm;

    // pktidx_stride is the size of a single channel for a single PKTIDX value
    // (i.e. for a single packet):
    const int pktidx_stride = ata_oi->pkt_nchan;

    // Stream is the "channel chunk" for this FID
    const int stream = (feng_id - ata_oi->schan) / ata_oi->pkt_nchan;

    // Advance dst_base to...
    const long offset = feng_id * fid_stride // first location of this FID, then
            +  stream * stream_stride // first location of this stream, then
            +  (pkt_idx - d->pktidx_per_block) * pktidx_stride; // to this pktidx

    // if (! copy_packet_data_to_databuf_printed){
    //     printf("nstrm            = %d\n", ata_oi->nstrm);
    //     printf("pkt_nchan        = %d\n", ata_oi->pkt_nchan);
    //     printf("feng_id          = %d\n", feng_id);
    //     printf("schan            = %d\n", ata_oi->schan);
    //     printf("pkt_payload_size = %d\n", pkt_payload_size);
    //     printf("pktidx_per_block = %d\n", d->pktidx_per_block);
    //     printf("stream_stride    = %d\n", stream_stride);
    //     printf("fid_stride       = %d\n", fid_stride);
    //     printf("pktidx_stride    = %d\n", pktidx_stride);
    //     printf("stream           = %d\n", stream);
    //     printf("offset           = %ld\n", offset);
    //     copy_packet_data_to_databuf_printed = 1;
    // }

    dst_base += offset;


    /* Packet has full data, just do a memcpy */
    memcpy(dst_base, ata_pkt->payload, pkt_payload_size);
}

/* Push all blocks down a level, losing the first one */
static void block_stack_push(struct datablock_stats *d, int nblock)
{
    int i;
    for (i=1; i<nblock; i++)
        memcpy(&d[i-1], &d[i], sizeof(struct datablock_stats));
}

// Check the given pktidx value against the status buffer's PKTSTART/PKTSTOP
// values. Logic goes something like this:
//   if PKTSTART <= pktidx < PKTSTOPs
//     if STTVALID == 0
//       STTVALID=1
//       calculate and store STT_IMJD, STT_SMJD
//     endif
//     return RECORD
//   else
//     STTVALID=0
//     return ARMED
//   endif
static
enum run_states check_start_stop(hashpipe_status_t *st, uint64_t pktidx)
{
  enum run_states retval = ARMED;
  uint32_t sttvalid = 0;
  uint64_t pktstart = 0;
  uint64_t pktstop = 0;

  uint32_t pktntime = ATASNAP_DEFAULT_PKTNTIME;
  uint64_t synctime = 0;
  double chan_bw = 1.0;

  double realtime_secs = 0.0;
  struct timespec ts;

  int    stt_imjd = 0;
  int    stt_smjd = 0;
  double stt_offs = 0;

  hashpipe_status_lock_safe(st);
  {
    hgetu4(st->buf, "STTVALID", &sttvalid);
    hgetu8(st->buf, "PKTSTART", &pktstart);
    hgetu8(st->buf, "PKTSTOP", &pktstop);

    if(pktstart <= pktidx && pktidx < pktstop) {
      retval = RECORD;
      hputs(st->buf, "DAQSTATE", "RECORD");

      if(sttvalid != 1) {
        hputu4(st->buf, "STTVALID", 1);

        hgetu4(st->buf, "PKTNTIME", &pktntime);
        hgetr8(st->buf, "CHAN_BW", &chan_bw);
        hgetu8(st->buf, "SYNCTIME", &synctime);

        // Calc real-time seconds since SYNCTIME for pktidx:
        //
        //                      pktidx * pktntime
        //     realtime_secs = -------------------
        //                        1e6 * chan_bw
        if(chan_bw != 0.0) {
          realtime_secs = pktidx * pktntime / (1e6 * fabs(chan_bw));
        }

        ts.tv_sec = (time_t)(synctime + rint(realtime_secs));
        ts.tv_nsec = (long)((realtime_secs - rint(realtime_secs)) * 1e9);

        get_mjd_from_timespec(&ts, &stt_imjd, &stt_smjd, &stt_offs);

        hputu4(st->buf, "STT_IMJD", stt_imjd);
        hputu4(st->buf, "STT_SMJD", stt_smjd);
        hputr8(st->buf, "STT_OFFS", stt_offs);
      }
    } else {
      hputs(st->buf, "DAQSTATE", "ARMED");
      if(sttvalid != 0) {
        hputu4(st->buf, "STTVALID", 0);
      }
    }
  }
  hashpipe_status_unlock_safe(st);

  return retval;
}

static int init(hashpipe_thread_args_t *args)
{
  	const char * thread_name = args->thread_desc->name;

    /* Non-network essential paramaters */
    int blocsize=BLOCK_DATA_SIZE;
    int directio=1;
    int nbits=4;
    int npol=2;
    double obsbw=187.5;
    int obsnchan=64;
    int obsschan=0;
    int overlap=0;
    double tbin=0.0;
    char obs_mode[80] = {0};
    char fifo_name[PATH_MAX];

    /* Create control FIFO (/tmp/hpguppi_daq_control/$inst_id) */
    int rv = mkdir(HPGUPPI_DAQ_CONTROL, 0777);
    if (rv!=0 && errno!=EEXIST) {
        hashpipe_error(thread_name, "Error creating control fifo directory");
        return HASHPIPE_ERR_SYS;
    } else if(errno == EEXIST) {
        errno = 0;
    }

    sprintf(fifo_name, "%s/%d", HPGUPPI_DAQ_CONTROL, args->instance_id);
    rv = mkfifo(fifo_name, 0666);
    if (rv!=0 && errno!=EEXIST) {
        hashpipe_error(thread_name, "Error creating control fifo");
        return HASHPIPE_ERR_SYS;
    } else if(errno == EEXIST) {
        errno = 0;
    }

    struct hpguppi_pktsock_params *p_psp = (struct hpguppi_pktsock_params *)
        malloc(sizeof(struct hpguppi_pktsock_params));

    if(!p_psp) {
        perror(__FUNCTION__);
        return -1;
    }

    strcpy(obs_mode, "RAW");

    hashpipe_status_t *st = &args->st;

    hashpipe_status_lock_safe(st);
    // Get network parameters (BINDHOST, BINDPORT, PKTFMT)
    hpguppi_read_pktsock_params(st->buf, p_psp);
    // Get info from status buffer if present (no change if not present)
    hgeti4(st->buf, "BLOCSIZE", &blocsize);
    hgeti4(st->buf, "DIRECTIO", &directio);
    hgeti4(st->buf, "NBITS", &nbits);
    hgeti4(st->buf, "NPOL", &npol);
    hgetr8(st->buf, "OBSBW", &obsbw);
    hgeti4(st->buf, "OBSNCHAN", &obsnchan);
    hgeti4(st->buf, "OBSSCHAN", &obsschan);
    hgeti4(st->buf, "OVERLAP", &overlap);
    hgets(st->buf, "OBS_MODE", sizeof(obs_mode), obs_mode);
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
    hputr8(st->buf, "OBSBW", obsbw);
    hputi4(st->buf, "OBSNCHAN", obsnchan);
    hputi4(st->buf, "OBSSCHAN", obsschan);
    hputi4(st->buf, "OVERLAP", overlap);
    hputr8(st->buf, "TBIN", tbin);
    hputs(st->buf, "OBS_MODE", obs_mode);
    // Data are in time-major order (i.e. time dimension changes faster than
    // channel dimension), so specify that data are NOT in channel major order.
    hputi4(st->buf, "CHANMAJ", 0);
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

    rv = hashpipe_pktsock_open(&p_psp->ps, p_psp->ifname, PACKET_RX_RING);
    if (rv!=HASHPIPE_OK) {
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

    /* Open command FIFO for read */
    char fifo_name[PATH_MAX];
    char fifo_cmd[MAX_CMD_LEN];
    sprintf(fifo_name, "%s/%d", HPGUPPI_DAQ_CONTROL, args->instance_id);
    int fifo_fd = open(fifo_name, O_RDONLY | O_NONBLOCK);
    if (fifo_fd<0) {
        hashpipe_error(thread_name, "Error opening control fifo)");
        pthread_exit(NULL);
    }

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

    /* Time parameters */
    // int stt_imjd=0, stt_smjd=0;
    // double stt_offs=0.0;

    /* Packet format to use */
    int seq_step = 1;
    int pkt_nchan=0, pkt_npol=0, pkt_nbits=0, pkt_payload_size;
    pkt_nchan = pf.hdr.nchan;
    pkt_nbits = pf.hdr.nbits;
    pkt_npol = pf.hdr.npol;
    pkt_payload_size = pkt_nchan*16*pkt_npol*(pkt_nbits*2/8);//exclude the ATA-SNAP's header-16 bytes 
    size_t packet_data_size = pkt_payload_size+16;//hpguppi_udp_packet_datasize(p_ps_params->packet_size);

    fprintf(stderr, "Observational header:\n"
            "\tnchan: %d\n"
            "\tnbits: %d\n"
            "\tnpol: %d\n"
            "\t\tPacket Payload Size (nchan*16*npol*(nbits*2/8)): %d\n",
            pkt_nchan,
            pkt_nbits,
            pkt_npol,
            pkt_payload_size
    );

    // Capture expected step change in ATA SNAP packet's sequential timestamps
    hashpipe_status_lock_safe(st);
    hgeti4(status_buf, "PKTNTIME", &seq_step);
    hashpipe_status_unlock_safe(st);

    /* Figure out size of data in each packet, number of packets
     * per block, etc.  Changing packet size during an obs is not
     * recommended.
     */
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
    unsigned pkt_per_block = block_size / pkt_payload_size;
    unsigned pktidx_per_block = pkt_per_block/seq_step;//ata_snap_pktidx_per_block(BLOCK_DATA_SIZE, obs_info);
    fprintf(stderr, "Packets per block %d, Packet timestamps per block %d\n", pkt_per_block, pktidx_per_block);
    unsigned long pkt_blk_num;


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
    unsigned i;
    const int n_wblock = 2;
    struct datablock_stats wblk[n_wblock];
    // Initialize working blocks
    for(i=0; i<n_wblock; i++) {
        //important to initialise the blocks with unique block idx an block_num values
        init_datablock_stats(&wblk[i], db, i, i, pkt_per_block);
        wait_for_block_free(wblk+i, st, status_key);
    }

    /* Misc counters, etc */
    int rv;
    // char *curdata=NULL, *curheader=NULL;
    uint64_t start_pkt_seq_num=0, stop_pkt_seq_num=0, pkt_seq_num, last_pkt_seq_num=2048, nextblock_pkt_seq_num=0;
    int64_t pkt_seq_num_diff;    
    // double drop_frac_avg=0.0;
    // const double drop_lpf = 0.25;
    // int netbuf_full = 0;
    // char netbuf_status[128] = {};
    // unsigned force_new_block=0, waiting=-1;
    unsigned waiting=-1;
    enum run_states state = IDLE;

    // Heartbeat variables
    time_t lasttime = 0;
    time_t curtime = 0;
    char timestr[32] = {0};
    
    unsigned int packets_per_second = 0;
    unsigned int packets_per_second_buf[PKTS_PER_SECOND_BUF_LENGTH] = {0};
    float average_packets_per_second = 0.0;

    // Drop all packets to date
    unsigned char *p_frame;
    while((p_frame=hashpipe_pktsock_recv_frame_nonblock(&p_ps_params->ps))) {
        hashpipe_pktsock_release_frame(p_frame);
    }
    struct ata_snap_pkt *ata_snap_pkt;

    // Get any obs info from status buffer, store values
    hashpipe_status_lock_safe(st);
    {
        // Read (no change if not present)
        hgetu4(st->buf, "FENCHAN",  &obs_info.fenchan);
        hgetu4(st->buf, "NANTS",    &obs_info.nants);
        hgetu4(st->buf, "NSTRM",    &obs_info.nstrm);
        hgetu4(st->buf, "PKTNTIME", &obs_info.pkt_ntime);
        hgetu4(st->buf, "PKTNCHAN", &obs_info.pkt_nchan);
        hgeti4(st->buf, "SCHAN",    &obs_info.schan);
        // If obs_info is valid
        if(ata_snap_obs_info_valid(obs_info)) {
        // Update obsnchan, pktidx_per_block, and eff_block_size
        // obsnchan = ata_snap_obsnchan(obs_info);
        // pktidx_per_block = ata_snap_pktidx_per_block(BLOCK_DATA_SIZE, obs_info);
        // eff_block_size = ata_snap_block_size(BLOCK_DATA_SIZE, obs_info);

        hputs(st->buf, "OBSINFO", "VALID");
        } else {
        hputs(st->buf, "OBSINFO", "INVALID");
        }

        // Write (store default/invlid values if not present)
        hputu4(st->buf, "FENCHAN",  obs_info.fenchan);
        hputu4(st->buf, "NANTS",    obs_info.nants);
        hputu4(st->buf, "NSTRM",    obs_info.nstrm);
        hputu4(st->buf, "PKTNTIME", obs_info.pkt_ntime);
        hputu4(st->buf, "PKTNCHAN", obs_info.pkt_nchan);
        hputi4(st->buf, "SCHAN",    obs_info.schan);

        // hputu4(st->buf, "OBSNCHAN", obsnchan);
        hputu4(st->buf, "PIPERBLK", pktidx_per_block);
        // hputi4(st->buf, "BLOCSIZE", eff_block_size);
    }
    hashpipe_status_unlock_safe(st);

    fprintf(stderr, "Receiving at interface %s, port %d, expecting packet size %ld\n",
    p_ps_params->ifname, p_ps_params->port, packet_data_size);//p_ps_params->packet_size);

    /* Main loop */
    while (run_threads()) {

        /* Wait for data */
        do {
            p_frame = hashpipe_pktsock_recv_udp_frame(
                &p_ps_params->ps, p_ps_params->port, 1000); // 1 second timeout
            
            // Heartbeat update?
            time(&curtime);
            if(curtime != lasttime) {//time stores seconds since epoch
                lasttime = curtime;
                average_packets_per_second = 0.0;
                for (int pps = 1; pps < PKTS_PER_SECOND_BUF_LENGTH; pps++){
                    packets_per_second_buf[pps-1] = packets_per_second_buf[pps];
                    average_packets_per_second += packets_per_second_buf[pps];
                }
                packets_per_second_buf[PKTS_PER_SECOND_BUF_LENGTH-1] = packets_per_second;
                packets_per_second = 0;
                average_packets_per_second = (packets_per_second + average_packets_per_second)/PKTS_PER_SECOND_BUF_LENGTH;

                ctime_r(&curtime, timestr);
                timestr[strlen(timestr)-1] = '\0'; // Chop off trailing newline
                hashpipe_status_lock_safe(st);
                {
                    hgetu4(st->buf, "FENCHAN",  &obs_info.fenchan);
                    hgetu4(st->buf, "NANTS",    &obs_info.nants);
                    hgetu4(st->buf, "NSTRM",    &obs_info.nstrm);
                    hgetu4(st->buf, "PKTNTIME", &obs_info.pkt_ntime);
                    hgetu4(st->buf, "PKTNCHAN", &obs_info.pkt_nchan);
                    hgeti4(st->buf, "SCHAN",    &obs_info.schan);
                    // If obs_info is valid
                    if(ata_snap_obs_info_valid(obs_info)) {
                        // Update obsnchan, pktidx_per_block, and eff_block_size
                        // obsnchan = ata_snap_obsnchan(obs_info);
                        // pktidx_per_block = ata_snap_pktidx_per_block(BLOCK_DATA_SIZE, obs_info);

                        // hputu4(st->buf, "OBSNCHAN", obsnchan);
                        hputu4(st->buf, "PIPERBLK", pktidx_per_block);
                        // hputi4(st->buf, "BLOCSIZE", eff_block_size);

                        hputs(st->buf, "OBSINFO", "VALID");
                    } else {
                        hputs(st->buf, "OBSINFO", "INVALID");
                    }

                    hputi8(st->buf, "PKTIDX", pkt_seq_num);
                    hputi8(st->buf, "NPKTS", npacket_total);
                    hputr4(st->buf, "PHYSPKPS", average_packets_per_second);
                    hputr4(st->buf, "PHYSGBPS", average_packets_per_second*packet_data_size/1e9);
                    hputs(st->buf, "DAQPULSE", timestr);
                    hputs(st->buf, "DAQSTATE", state == IDLE  ? "idle" :
                                            state == ARMED ? "armed"   : "record");
                    hputi8(st->buf, "NXTBLKSQN", nextblock_pkt_seq_num);
                }
                hashpipe_status_unlock_safe(st);
            }

            /* Set "waiting" flag */
            if (!p_frame && run_threads() && waiting!=1) {
                hashpipe_status_lock_safe(st);
                hputs(st->buf, status_key, "waiting");
                hashpipe_status_unlock_safe(st);
                waiting=1;
            }

            // Check FIFO for command
            rv = read(fifo_fd, fifo_cmd, MAX_CMD_LEN-1);
            if(rv == -1 && errno != EAGAIN) {
                hashpipe_error(thread_name, "error reading control fifo)");
            } else if(rv > 0) {
                // Trim newline from command, if any
                char *newline = strchr(fifo_cmd, '\n');
                if (newline!=NULL) *newline='\0';

                // Log command
                hashpipe_warn(thread_name, "got %s command", fifo_cmd);

                // Act on command
                if(strcasecmp(fifo_cmd, "QUIT") == 0) {
                    // Go to IDLE state
                    state = IDLE;
                    // Hashpipe will exit upon thread exit
                    pthread_exit(NULL);
                } else if(strcasecmp(fifo_cmd, "MONITOR") == 0) {
                    hashpipe_warn(thread_name,
                            "MONITOR command not supported, use null_output_thread.");
                } else if(strcasecmp(fifo_cmd, "START") == 0) {
                    // If in the IDLE or ARMED states
                    // (START in RECORD state is a no-op)
                    if(state == IDLE || state == ARMED) {
                        // Reset current blocks' packet_idx values and stats
                        for (i=0; i<n_wblock; i++) {
                            wblk[i].packet_idx = 0;
                            reset_datablock_stats(&wblk[i]);
                        }

                        // Go to (or stay in) ARMED state
                        state = ARMED;
                    }
                } else if(strcasecmp(fifo_cmd, "STOP") == 0) {
                    // Go to IDLE state
                    state = IDLE;
                } else {
                    hashpipe_error(thread_name,
                            "got unrecognized command '%s'", fifo_cmd);
                }
            }

        } while (!p_frame && run_threads());

        if(!run_threads()) {
            // We're outta here!
            hashpipe_pktsock_release_frame(p_frame);
            break;
        }

        // fprintf(stderr, "Packet received with size %d expecting size %ld\r", PKT_UDP_SIZE(p_frame) - 8, p_ps_params->packet_size);
        /* Check packet size */
        // if(p_ps_params->packet_size == 0) {
        //     p_ps_params->packet_size = PKT_UDP_SIZE(p_frame) - 8;
        // } else 
        // if(p_ps_params->packet_size != PKT_UDP_SIZE(p_frame) - 8) {
        if(packet_data_size != PKT_UDP_SIZE(p_frame) - 8) {
            /* Unexpected packet size, ignore? */
            nbogus_total++;
            if(nbogus_total % 10 == 0) {
                hashpipe_status_lock_safe(st);
                hputi4(st->buf, "NBOGUS", nbogus_total);
                hputi4(st->buf, "PKTSIZE", PKT_UDP_SIZE(p_frame)-8);
                hashpipe_status_unlock_safe(st);
            }
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
        
        packets_per_second++;

        ata_snap_pkt = (struct ata_snap_pkt*) p_frame;
        // Get packet's sequence number
        pkt_seq_num = __bswap_64(ata_snap_pkt->timestamp);

        // Update PKTIDX in status buffer if pkt_seq_num % pktidx_per_block == 0
        // and read PKTSTART, DWELL to calculate start/stop seq numbers.
        if(pkt_seq_num % pktidx_per_block == 0 || pkt_seq_num >= stop_pkt_seq_num) {// || pkt_seq_num < stop_pkt_seq_num - pktidx_per_block) {
            start_pkt_seq_num = pkt_seq_num;

            hashpipe_status_lock_safe(st);
            hputi8(st->buf, "PKTIDX", pkt_seq_num);
            // hgetu8(st->buf, "PKTSTART", &start_pkt_seq_num);
            // start_pkt_seq_num -= start_pkt_seq_num % pktidx_per_block;
            hputu8(st->buf, "PKTSTART", start_pkt_seq_num);
            // hgetr8(st->buf, "DWELL", &dwell_seconds);
            // // Dwell blocks is equal to:
            // //
            // //           dwell_seconds * samples/second
            // //     -------------------------------------------
            // //     samples/spectrum * spectra/pkt * pkts/block
            // //
            // // To get an integer number of blocks, simply truncate
            // dwell_blocks = trunc(dwell_seconds * samples_per_second /
            //         (samples_per_spectrum * spectra_per_packet * pktidx_per_block));

            stop_pkt_seq_num = start_pkt_seq_num + pktidx_per_block;// * dwell_blocks;
            hputi8(st->buf, "PKTSTOP", stop_pkt_seq_num);
            hashpipe_status_unlock_safe(st);
        }

        // If IDLE, release frame and continue main loop
        // if(state == IDLE) {
        //     // Release frame!
        //     hashpipe_pktsock_release_frame(p_frame);
        //     continue;
        // }

        /* Check seq num diff */
        pkt_seq_num_diff = pkt_seq_num - last_pkt_seq_num;
        if (pkt_seq_num_diff<seq_step) {
            if (pkt_seq_num_diff<-1024) {
                // force_new_block=1;
            } else if (pkt_seq_num_diff==0) {
                hashpipe_warn(thread_name,
                        "Received duplicate packet (pkt_seq_num=%lu)",
                        pkt_seq_num);
            }
            else {
              // Release frame!
              hashpipe_pktsock_release_frame(p_frame);
              /* No going backwards */
              continue;
            }
        } else {
            // force_new_block=0;
            if(npacket_total == 0) {
                npacket_total = 1;
                ndrop_total = 0;
            } else {
                npacket_total += pkt_seq_num_diff/seq_step;
                ndrop_total += (pkt_seq_num_diff - seq_step)/seq_step;
            }
        }

        // // Ignore packets with FID >= NANTS
        // if(ata_snap_pkt->feng_id >= obs_info.nants) {
        //     continue;
        // }

        pkt_blk_num = pkt_seq_num / pktidx_per_block;
        // Manage blocks based on pkt_blk_num
        if(pkt_blk_num == wblk[n_wblock-1].block_num + 1) {
            // Time to advance the blocks!!!
            // Finalize first working block
            finalize_block(wblk);
            // Update ndrop counter
            ndrop_total += wblk->ndrop;
            // Shift working blocks
            block_stack_push(wblk, n_wblock);
            // Check start/stop using wblk[0]'s first PKTIDX
            check_start_stop(st, pkt_seq_num);//wblk[0].block_num * pktidx_per_block);
            // Increment last working block
            increment_block(&wblk[n_wblock-1], pkt_blk_num);
            // Wait for new databuf data block to be free
            wait_for_block_free(&wblk[n_wblock-1], st, status_key);
        }
        // Check for PKTIDX discontinuity
        else if(pkt_blk_num < wblk[0].block_num - 1
                || pkt_blk_num > wblk[n_wblock-1].block_num + 1) {
            // Should only happen when transitioning into ARMED, so warn about it
            hashpipe_warn(thread_name,
                "working blocks reinit due to packet discontinuity (PKTIDX %lu)",
                pkt_seq_num);

            // Re-init working blocks for block number *after* current packet's block
            // and clear their data buffers
            for(i=0; i<n_wblock; i++) {
            init_datablock_stats(wblk+i, NULL, -1, pkt_blk_num+i+1,
                pkt_per_block);
            }
            fprintf(stderr, "Packet Block Index for finalisation: %ld\n", wblk[n_wblock-1].block_num + 1);
            // Check start/stop using wblk[0]'s first PKTIDX
            check_start_stop(st, wblk[0].block_num * pktidx_per_block);
            // This happens after discontinuities (e.g. on startup), so don't warn about
            // it.
        } else if(pkt_blk_num == wblk[0].block_num - 1) {
            // Ignore late packet, continue on to next one
            // TODO Move this check above the "once per block" status buffer
            // update (so we don't accidentally update status buffer based on a
            // late packet)?
            // nlate++;
        }

        // TODO Check START/STOP status???

        // Once we get here, compute the index of the working block corresponding
        // to this packet.  The computed index may not correspond to a valid
        // working block!
        int wblk_idx = pkt_blk_num - wblk[0].block_num;

        // Only copy packet data and count packet if its wblk_idx is valid
        if(0 <= wblk_idx && wblk_idx < n_wblock) {
            // Update block's packets per block.  Not needed for each packet, but
            // probably just as fast to do it for each packet rather than
            // check-and-update-only-if-needed for each packet.
            wblk[wblk_idx].pkts_per_block = pkt_per_block;//eff_block_size / ATA_SNAP_PKT_SIZE_PAYLOAD;
            wblk[wblk_idx].pktidx_per_block = pktidx_per_block;

            // Copy packet data to data buffer of working block
            copy_packet_data_to_databuf(wblk+wblk_idx,
                &obs_info, ata_snap_pkt, pkt_payload_size);

            // Count packet for block and for processing stats
            wblk[wblk_idx].npacket++;
        }


        last_pkt_seq_num = pkt_seq_num;
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

// vi: set ts=8 sw=4 et :
