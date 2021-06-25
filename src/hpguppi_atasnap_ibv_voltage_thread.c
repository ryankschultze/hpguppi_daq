// hpguppi_atasnap_ibv_voltage_thread.c
//
// A Hashpipe thread that proceeses "voltage mode" packets sent by the ATA SNAP
// design from an input buffer (populated by hpguppi_ibverbs_pkt_thread) and
// assembles them into GUPPI RAW blocks.

// TODO TEST Wait for first (second?) start-of-block when transitioning into
//           LISTEN state so that the first block will be complete.
// TODO Add PSPKTS and PSDRPS status buffer fields for pktsock
// TODO TEST Set NETSTAE to idle in IDLE state
// TODO TEST IP_DROP_MEMBERSHIP needs mcast IP address (i.e. not 0.0.0.0)

#define _GNU_SOURCE 1
//#include <stdio.h>
//#include <sys/types.h>
#include <stdlib.h>
#include <sched.h>
#include <math.h>
#include <unistd.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <poll.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/time.h>
#include <sys/resource.h>

#include "hashpipe.h"
#include "hpguppi_databuf.h"
#include "hpguppi_time.h"
#include "hpguppi_util.h"
#include "hpguppi_atasnap.h"
#include "hpguppi_ibverbs_pkt_thread.h"

// Change to 1 to use temporal memset() rather than non-temporal bzero_nt()
#if 0
#define bzero_nt(d,l) memset(d,0,l)
#endif

// Change to 1 to use temporal memcpy() rather than non-temporal memcpy_nt()
#if 0
#define memcpy_nt(dst,src,len) memcpy(dst,src,len)
#endif

#define ELAPSED_NS(start,stop) \
  (((int64_t)stop.tv_sec-start.tv_sec)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

#define HPUT_DAQ_STATE(st, state)\
  hputs(st->buf, "DAQSTATE", state == IDLE  ? "idling" :\
                             state == ARMED ? "armed"  :\
                             "recording")

// This thread's init() function, if provided, is called by the Hashpipe
// framework at startup to allow the thread to perform initialization tasks
// such as setting up network connections or GPU devices.
static int init(hashpipe_thread_args_t *args)
{
  // Local aliases to shorten access to args fields
  // Our input buffer happens to be a hpguppi_input_databuf
  hpguppi_input_databuf_t *dbin  = (hpguppi_input_databuf_t *)args->ibuf;
  const char * thread_name = args->thread_desc->name;
  const char * status_key = args->thread_desc->skey;
  hashpipe_status_t *st = &args->st;

  // Non-network essential paramaters
  int blocsize=BLOCK_DATA_SIZE;
  int directio=1;
  int nbits=4;
  int npol=4;
  double obsfreq=0;
  double chan_bw=900.0/4096;
  double obsbw=256*chan_bw;
  int obsnchan=1;
  int nants=1;
  int overlap=0;
  double tbin=1e-6;
  char obs_mode[80] = {0};
  struct rlimit rlim;

  strcpy(obs_mode, "RAW");

  // Verify that the IBVPKTSZ was specified as expected/requried
  if(hpguppi_pktbuf_slot_offset(dbin, ATA_SNAP_PKT_OFFSET_HEADER) %
      PKT_ALIGNMENT_SIZE != 0
  || hpguppi_pktbuf_slot_offset(dbin, ATA_SNAP_PKT_OFFSET_PAYLOAD) %
      PKT_ALIGNMENT_SIZE != 0) {
    errno = EINVAL;
    hashpipe_error(thread_name, "IBVPKTSZ!=%d,%d,[...]",
        ATA_SNAP_PKT_OFFSET_HEADER, ATA_SNAP_PKT_SIZE_HEADER);
    return HASHPIPE_ERR_PARAM;
  }

  // Set RLIMIT_RTPRIO to 1
  getrlimit(RLIMIT_RTPRIO, &rlim);
	rlim.rlim_cur = 1;
	rlim.rlim_max = (rlim.rlim_max < rlim.rlim_cur ? rlim.rlim_cur : rlim.rlim_max);
  if(setrlimit(RLIMIT_RTPRIO, &rlim)) {
    hashpipe_error(thread_name, "setrlimit(RLIMIT_RTPRIO)");
    return HASHPIPE_ERR_PARAM;
  }

  struct sched_param sched_param = {
    .sched_priority = 1
  };
  if(sched_setscheduler(0, SCHED_RR, &sched_param)) {
    hashpipe_error(thread_name, "sched_setscheduler");
  }

  hashpipe_status_lock_safe(st);
  {
    // Get info from status buffer if present (no change if not present)
    hgeti4(st->buf, "BLOCSIZE", &blocsize);
    hgeti4(st->buf, "DIRECTIO", &directio);
    hgeti4(st->buf, "NANTS", &nants);
    hgeti4(st->buf, "NBITS", &nbits);
    hgeti4(st->buf, "NPOL", &npol);
    hgetr8(st->buf, "OBSFREQ", &obsfreq);
    hgetr8(st->buf, "OBSBW", &obsbw);
    hgetr8(st->buf, "CHAN_BW", &chan_bw);
    hgeti4(st->buf, "OBSNCHAN", &obsnchan);
    hgeti4(st->buf, "OVERLAP", &overlap);
    hgets(st->buf, "OBS_MODE", sizeof(obs_mode), obs_mode);

    // Prevent div-by-zero errors (should never happen...)
    if(nants == 0) {
      nants = 1;
      hputi4(st->buf, "NANTS", nants);
    }

    // If CHAN_BW is zero, set to default value (1 MHz)
    if(chan_bw == 0.0) {
      chan_bw = 1.0;
    }

    // Calculate tbin and obsbw from chan_bw
    tbin = 1e-6 / fabs(chan_bw);
    obsbw = chan_bw * obsnchan / nants;

    // Update status buffer (in case fields were not there before).
    hputs(st->buf, "DAQSTATE", "LISTEN");
    hputi4(st->buf, "BLOCSIZE", blocsize);
    hputi4(st->buf, "DIRECTIO", directio);
    hputi4(st->buf, "NBITS", nbits);
    hputi4(st->buf, "NPOL", npol);
    hputr8(st->buf, "OBSBW", obsbw);
    hputr8(st->buf, "CHAN_BW", chan_bw);
    hputi4(st->buf, "OBSNCHAN", obsnchan);
    hputi4(st->buf, "OVERLAP", overlap);
    hputs(st->buf, "PKTFMT", "ATASNAPV");
    hputr8(st->buf, "TBIN", tbin);
    hputs(st->buf, "OBS_MODE", obs_mode);
    hputi4(st->buf, "NDROP", 0);
    // Set status_key to init
    hputs(st->buf, status_key, "init");
  }
  hashpipe_status_unlock_safe(st);

  // Success!
  return 0;
}

static void * run(hashpipe_thread_args_t * args)
{
#if 0
int debug_i=0, debug_j=0;
#endif
  // Local aliases to shorten access to args fields
  // Our input and output buffers happen to be a hpguppi_input_databuf
  hpguppi_input_databuf_t *dbin  = (hpguppi_input_databuf_t *)args->ibuf;
  hpguppi_input_databuf_t *dbout = (hpguppi_input_databuf_t *)args->obuf;
  hashpipe_status_t *st = &args->st;
  const char * thread_name = args->thread_desc->name;
  const char * status_key = args->thread_desc->skey;

  // Max flows allowed (optionally from hpguppi_ibvpkt_thread via status
  // buffer)
  uint32_t max_flows = 16;
  // Port to listen on
  uint32_t port = 10000;

  // Current run state
  //enum run_states state = LISTEN;
  unsigned waiting = 0;
  // Update status_key with idle state and get max_flows, port
  hashpipe_status_lock_safe(st);
  {
    hputs(st->buf, status_key, "listen");
    hgetu4(st->buf, "MAXFLOWS", &max_flows);
    hgetu4(st->buf, "BINDPORT", &port);
    // Store bind port in status buffer (in case it was not there before).
    hputu4(st->buf, "BINDPORT", port);
  }
  hashpipe_status_unlock_safe(st);

  // Make sure we got a non-zero max_flows
  if(max_flows == 0) {
    hashpipe_error(thread_name, "MAXFLOWS must be non-zero!");
    return NULL;
  }

  // Misc counters, etc
  int rv=0;
  int i;

#if 0
  uint64_t u64;
  uint8_t u8 = 0;
  uint8_t *pu8in = (uint8_t *)dbin;
  uint8_t *pu8out = (uint8_t *)dbout;
  for(u64=0; u64<sizeof(hpguppi_input_databuf_t); u64+=4096) {
    if(u8 || !u8) {
      u8 += pu8in[u64];
      u8 += pu8out[u64];
    }
  }
  hashpipe_info(thread_name, "db pagein sum is %u", u8);
#endif
  memset(dbout->block, 0, sizeof(dbout->block));
  hashpipe_info(thread_name,
      "set %lu bytes in dbout to 0", sizeof(dbout->block));

  // for(i=0; i<N_INPUT_BLOCKS; i++) {
  //   hashpipe_info(thread_name, "db_in  block %2d : %p %p", i,
  //       hpguppi_databuf_data(dbin, i),
  //       hpguppi_databuf_data(dbin, i) + BLOCK_DATA_SIZE - 1);
  // }

  // for(i=0; i<N_INPUT_BLOCKS; i++) {
  //   hashpipe_info(thread_name, "db_out block %2d : %p %p", i,
  //       hpguppi_databuf_data(dbout, i),
  //       hpguppi_databuf_data(dbout, i) + BLOCK_DATA_SIZE - 1);
  // }

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
  // wblk is a two element array of block_info structures (i.e. the working
  // blocks)
  int wblk_idx;
  const int n_wblock = 2;
  struct datablock_stats wblk[n_wblock];

  // Packet block variables
  uint64_t pkt_seq_num = 0;
  unsigned long pkt_blk_num, last_pkt_blk_num = ~0;
  uint64_t obs_start_seq_num = 0, obs_stop_seq_num = 0, pkt_obs_relative_idx = 0;
  uint64_t status_seq_num;
  uint64_t pkt_payload_size;
  uint32_t fid_stride, time_stride;
  uint32_t pkt_stream;

  // Heartbeat variables
  struct timespec ts_start_block = {0}, ts_stop_block = {0};
  time_t lasttime = 0;
  time_t curtime = 0;
  char timestr[32] = {0};
  float blocks_per_second = 0.0;

  // State variables
  enum run_states state = IDLE;
  char flag_state_update = 0;
  char  PKT_OBS_IDX_flagged, PKT_OBS_FENG_flagged, 
        PKT_OBS_SCHAN_flagged, PKT_OBS_STREAM_flagged, pkt_obs_code;
  char flag_obs_end = 0, flag_arm_for_obs = 0;

  // Variables for working with the input databuf
  struct hpguppi_pktbuf_info * pktbuf_info = hpguppi_pktbuf_info_ptr(dbin);
  int block_idx_in = 0;
  struct timespec timeout_in = {0, 50 * 1000 * 1000}; // 50 ms
  char * datablock_header;

  // Variables for counting packets and bytes.
  uint64_t npacket_total=0, ndrop_total=0, /* nbogus_total=0, */ lastnpacket_total=0;
  uint64_t obs_npacket_total=0, obs_ndrop_total=0, obs_block_discontinuities=0;
  // uint32_t nbogus_size;

  // Variables for handing received packets
  uint8_t * p_u8pkt;
  struct ata_snap_ibv_pkt * p_pkt = NULL;
  const uint8_t * p_payload = NULL;

  struct mjd_t *mjd = malloc(sizeof(struct mjd_t));

  // Structure to hold observation info, init all fields to invalid values
  struct ata_snap_obs_info obs_info;
  ata_snap_obs_info_init(&obs_info);

  // Structure to hold feng info from packet
  struct ata_snap_feng_info feng_info = {0};

  // Variables for tracking timing stats
  //
  // ts_start_recv(N) to ts_stop_recv(N) is the time spent in the "receive" call.
  // ts_stop_recv(N) to ts_start_recv(N+1) is the time spent processing received data.
  struct timespec ts_start_recv = {0}, ts_stop_recv = {0};
  struct timespec ts_input_full0 = {0};
  struct timespec ts_free_input = {0};

  // Used to calculate moving average of fill-to-free times for input blocks
  uint64_t fill_to_free_elapsed_ns;
  uint64_t fill_to_free_moving_sum_ns = 0;
  uint64_t fill_to_free_block_ns[N_INPUT_BLOCKS] = {0};

  //struct timespec ts_sleep = {0, 10 * 1000 * 1000}; // 10 ms

#if 0
  // Allocate a 2K buffer into which packet will be non-temporally copied
  // before processing.  This buffer will be cached (due to parsing of the
  // headers), but the input databuf blocks will not be cached.
  if((rv = posix_memalign((void **)&p_spdpkt, 4096, MAX_PKT_SIZE))) {
    errno = rv;
    hashpipe_error(thread_name, "cannot allocate page aligned packet buffer");
    return NULL;
  }
#endif

  // Initialize working blocks
  for(wblk_idx=0; wblk_idx<n_wblock; wblk_idx++) {
      // important to initialise the blocks with unique block idx
      // and to set the block_num < -1 so the first observation packet causes a 
      // discontinuity and so a contemporary re-init of the working blocks
      init_datablock_stats(wblk+wblk_idx, dbout, wblk_idx, -10, obs_info.pkt_per_block);
      wait_for_block_free(wblk+wblk_idx, st, status_key);
  }

  // Get any obs info from status buffer, store values
  ata_snap_obs_info_read(st, &obs_info);
  ata_snap_obs_info_write(st, &obs_info);

  // Wait for ibvpkt thread to be running, then it's OK to add/remove flows.
  hpguppi_ibvpkt_wait_running(st);

  uint8_t* mac = ((struct hashpipe_ibv_context *)(dbin->padding + sizeof(struct hpguppi_pktbuf_info)))->mac;
  hashpipe_info(thread_name, "DEST MAC: %02x:%02x:%02x:%02x:%02x:%02x", mac[0], mac[1], mac[2], mac[3], mac[4], mac[5]);
  if(hpguppi_ibvpkt_flow(dbin, 0, IBV_FLOW_SPEC_UDP,
        mac, NULL, 0, 0,
        0, 0, 0, port))
  {
    hashpipe_error(thread_name, "hashpipe_ibv_flow error: (errno %d)", errno);
    hashpipe_info(thread_name, "exiting!");
    pthread_exit(NULL);

    return NULL;
  }

  // Main loop
  while (run_threads()) {

    // Mark ts_stop_recv as unset
    ts_stop_recv.tv_sec = 0;

    // Wait for data
    do {
      clock_gettime(CLOCK_MONOTONIC, &ts_start_recv);
      // If ts_stop_recv has been set
      // if(ts_stop_recv.tv_sec != 0) {
      //   // Accumulate processing time
      //   ns_processed_net += ELAPSED_NS(ts_stop_recv, ts_start_recv);
      // }
      rv = hpguppi_input_databuf_wait_filled_timeout(
          dbin, block_idx_in, &timeout_in);
      clock_gettime(CLOCK_MONOTONIC, &ts_stop_recv);

      time(&curtime);

      if(rv == HASHPIPE_TIMEOUT && curtime == lasttime) {
        // No, continue receiving
        continue;
      }

      // Got packets or new second

      // We perform some status buffer updates every second
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
          // hputi4(st->buf, "NBOGUS", nbogus_total);
          // hputi4(st->buf, "BOGUSIZE", nbogus_size);
          hputs(st->buf, "DAQPULSE", timestr);
          HPUT_DAQ_STATE(st, state);
        }
        hashpipe_status_unlock_safe(st);
        lastnpacket_total = npacket_total;
      } // curtime != lasttime

      // Set status field to "waiting" if we are not getting packets
      if(rv && run_threads() && !waiting) {
        hashpipe_status_lock_safe(st);
        {
          hputs(st->buf, status_key, "waiting");
        }
        hashpipe_status_unlock_safe(st);
        waiting=1;
      }

      // Will exit if thread has been cancelled
      pthread_testcancel();
    } while (rv && run_threads()); // end wait for data loop

    if(!run_threads()) {
      // We're outta here!
      // But first mark the block free if we got one.
      if(!rv) {
        hpguppi_input_databuf_set_free(dbin, block_idx_in);
        clock_gettime(CLOCK_MONOTONIC, &ts_free_input);
        fprintf(stderr, "final fill-to-free %ld ns\n", ELAPSED_NS(ts_stop_recv, ts_free_input));
      }
      break;
    }

    //hashpipe_info(thread_name, "Got packets!");
    // hpguppi_input_databuf_set_free(dbin, block_idx_in);
    // block_idx_in = (block_idx_in + 1) % dbin->header.n_block;
    // continue;

    // Got packet(s)!  Update status if needed.
    if (waiting) {
      hashpipe_status_lock_safe(st);
      {
        hputs(st->buf, status_key, "receiving");
      }
      hashpipe_status_unlock_safe(st);
      waiting=0;
    }

    if(ts_input_full0.tv_sec == 0) {
      ts_input_full0 = ts_stop_recv;
    }

    // For each packet: process all packets
    p_u8pkt = (uint8_t *)hpguppi_databuf_data(dbin, block_idx_in);
    for(i=0; i < pktbuf_info->slots_per_block; i++, p_u8pkt += pktbuf_info->slot_size) {
#if 0
      // Non-temporally copy packet into cached buffer
      memcpy_nt(p_spdpkt, p_u8pkt, MAX_PKT_SIZE);
#else
      p_pkt = (struct ata_snap_ibv_pkt *)p_u8pkt;
#endif

      // Parse packet
      p_payload = ata_snap_parse_ibv_packet(p_pkt, &feng_info);

      /* Check packet size */
      // PKT_UDP_SIZE doesn't work on p_pkt, also the observability will be good enough...
      // if(obs_info.pkt_data_size != PKT_UDP_SIZE((unsigned char *)p_pkt) - 8) {
      //     /* Unexpected packet size, ignore */
      //     nbogus_total++;
      //     nbogus_size = PKT_UDP_SIZE((unsigned char *)p_pkt)-8;
      //     continue;
      // }

      // We update the status buffer at the start of each block
      // Also read PKTSTART and PKTSTOP
      if(pkt_seq_num % obs_info.pktidx_per_block == 0
          && pkt_seq_num != status_seq_num) {
        status_seq_num  = pkt_seq_num;

        hashpipe_status_lock_safe(st);
        {
          hputi8(st->buf, "PKTIDX", pkt_seq_num);

          hgetu8(st->buf, "PKTSTART", &obs_start_seq_num);
          // obs_start_seq_num -= obs_start_seq_num % obs_info.pktidx_per_block;
          // hputu8(st->buf, "PKTSTART", obs_start_seq_num);
          hgetu8(st->buf, "PKTSTOP", &obs_stop_seq_num);
        }
        hashpipe_status_unlock_safe(st);
      } // End status buffer block update

      // Get packet index and absolute block number for packet
      pkt_seq_num = feng_info.pktidx;
      pkt_stream = (feng_info.feng_chan - obs_info.schan) / obs_info.pkt_nchan;

      // if(pkt_obs_code != PKT_OBS_OK){
      //   hashpipe_info(thread_name, "pkt version: %u\n", p_pkt->version);
      //   hashpipe_info(thread_name, "pkt type: %u\n", p_pkt->type);
      //   hashpipe_info(thread_name, "pkt n_chans: %u\n", __bswap_16(p_pkt->n_chans));
      //   hashpipe_info(thread_name, "pkt chan: %u\n", feng_info.feng_chan);
      //   hashpipe_info(thread_name, "pkt feng_id: %u\n", feng_info.feng_id);
      //   hashpipe_info(thread_name, "pkt timestamp: %lu\n", feng_info.pktidx);
      // }

      switch(state_from_start_stop(pkt_seq_num, obs_start_seq_num, obs_stop_seq_num)){
          case IDLE:// If should IDLE, 
            if(state == ARMED){
              flag_state_update = 1;
              state = IDLE;
            }
            // if RECORDING -> IDLE handled before finalising a block
            break;
          case RECORD:// If should RECORD
            if (state != RECORD && ata_snap_obs_info_valid(obs_info)){// Only enter recording mode if obs_params are valid
              // and not recording, flag obs_start
              flag_state_update = 1;
              flag_arm_for_obs = (state == IDLE ? 1 : 0);
              state = RECORD;
            }
            break;
          case ARMED:// If should ARM,
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
          
          ata_snap_obs_info_read(st, &obs_info);
          obs_npacket_total = 0;
          obs_ndrop_total = 0;
          obs_block_discontinuities = 0;
          update_stt_status_keys(st, state, obs_start_seq_num, mjd);

          pkt_payload_size = ata_snap_pkt_payload_bytes(obs_info);
          fid_stride = obs_info.nstrm*pkt_payload_size;
          time_stride = obs_info.nants*fid_stride;

          for(wblk_idx=0; wblk_idx<n_wblock; wblk_idx++) {
            wblk[wblk_idx].pkts_per_block = obs_info.pkt_per_block;
            wblk[wblk_idx].pktidx_per_block = obs_info.pktidx_per_block;
            wblk[wblk_idx].packet_idx = obs_start_seq_num + wblk_idx*obs_info.pktidx_per_block;
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
          }
          hashpipe_info(thread_name, "Armed wblk for observation: obs_start_seq_num = %lu", obs_start_seq_num);
        }

        pkt_obs_relative_idx = pkt_seq_num - obs_start_seq_num;
        pkt_blk_num = pkt_obs_relative_idx / obs_info.pktidx_per_block;
        // fprintf(stderr, "seq: %012ld\tblk: %06ld\n", pkt_seq_num, pkt_blk_num);

        // Update PKTIDX in status buffer if pkt_seq_num % obs_info.pktidx_per_block == 0
        if(pkt_blk_num != last_pkt_blk_num){

          last_pkt_blk_num = pkt_blk_num;

          hashpipe_status_lock_safe(st);
          hputu8(st->buf, "BLKIDX", pkt_blk_num);
          hputu8(st->buf, "PKTIDX", pkt_seq_num);
          hashpipe_status_unlock_safe(st);
        }

        // Tally npacket_totals
        npacket_total += 1;
        // count ndrop only when a block is finalised, lest it is filled out of order.
        
        if (state != RECORD){
          continue;
        }
        
        pkt_obs_code = check_pkt_observability(&obs_info,// ata_snap_pkt, 
                    pkt_seq_num, obs_start_seq_num, feng_info.feng_id, pkt_stream, feng_info.feng_chan);

        // Manage blocks based on pkt_blk_num
        // fprintf(stderr, "%010ld\r", pkt_blk_num);
        if(pkt_blk_num == wblk[n_wblock-1].block_num + 1) {
            // Time to advance the blocks!!!
            // fprintf(stderr, "\nFinalising Block: %ld", wblk[0].block_num);
            clock_gettime(CLOCK_MONOTONIC, &ts_stop_block);

            // If the block to be finalised is out of observation range then flag_obs_end
            if(state_from_start_stop(wblk[0].packet_idx, obs_start_seq_num, obs_stop_seq_num) == IDLE &&
              state_from_start_stop(wblk[1].packet_idx, obs_start_seq_num, obs_stop_seq_num) == IDLE){
                  flag_obs_end = 1;
            }

            // Finalize first working block
            datablock_header = datablock_stats_header(&wblk[0]);
            // update PKTIDX triggering rawdisk to close fd when this last block gets finalised.
            hputu8(datablock_header, "PKTIDX", flag_obs_end ? pkt_seq_num : wblk[0].packet_idx);
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
              wblk[wblk_idx].packet_idx = pkt_seq_num + wblk_idx*obs_info.pktidx_per_block;
              init_datablock_stats(wblk+wblk_idx, NULL, -1,
                  pkt_blk_num+wblk_idx,
                  obs_info.pkt_per_block);

              // also update the working blocks' headers
              datablock_header = datablock_stats_header(wblk+wblk_idx);
              hashpipe_status_lock_safe(st);
                memcpy(datablock_header, st->buf, HASHPIPE_STATUS_TOTAL_SIZE);
              hashpipe_status_unlock_safe(st);
            }
            // Immediately update STT keys to ensure that the rawdisk thread uses a new stem
            update_stt_status_keys(st, state, obs_start_seq_num, mjd);
        }

        // Check observation state
        if(flag_obs_end){// transitioning RECORD->IDLE
          flag_obs_end = 0;
          obs_start_seq_num = 0;
          state = IDLE;
          flag_state_update = 1;
          update_stt_status_keys(st, state, pkt_seq_num, mjd);
          // Skip the reset of the incoming block
          break;
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
                  p_payload, pkt_obs_relative_idx,
                  feng_info.feng_id, pkt_stream, feng_info.feng_chan,
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

      if(pkt_obs_code != PKT_OBS_OK){
        // Skip the whole block of packets
        state = IDLE;
        flag_state_update = 1;
        break;
      }
    } // end for each packet

    // Mark input block free
    hpguppi_input_databuf_set_free(dbin, block_idx_in);

    // Update moving sum (for moving average)
    clock_gettime(CLOCK_MONOTONIC, &ts_free_input);
    fill_to_free_elapsed_ns = ELAPSED_NS(ts_stop_recv, ts_free_input);
    // Add new value, subtract old value
    fill_to_free_moving_sum_ns +=
        fill_to_free_elapsed_ns - fill_to_free_block_ns[block_idx_in];
    // Store new value
    fill_to_free_block_ns[block_idx_in] = fill_to_free_elapsed_ns;

    if(block_idx_in == N_INPUT_BLOCKS - 1) {
      hashpipe_status_lock_safe(st);
      {
        hputr8(st->buf, "NETBLKMS",
            round((double)fill_to_free_moving_sum_ns / N_INPUT_BLOCKS) / 1e6);
      }
      hashpipe_status_unlock_safe(st);
    }

#if 0
    fprintf(stderr, "blkin %d fill at %ld free +%ld ns (%d packets)\n",
        block_idx_in,
        ELAPSED_NS(ts_input_full0, ts_stop_recv),
        ELAPSED_NS(ts_stop_recv, ts_free_input), njobs);
#endif

    // Advance to next input block
    block_idx_in = (block_idx_in + 1) % dbin->header.n_block;

    // Will exit if thread has been cancelled
    pthread_testcancel();
  } // end main loop

  hashpipe_info(thread_name, "exiting!");
  pthread_exit(NULL);

  return NULL;
}

static hashpipe_thread_desc_t thread_desc = {
    name: "hpguppi_atasnap_ibv_voltage_thread",
    skey: "NETSTAT",
    init: init,
    run:  run,
    ibuf_desc: {hpguppi_input_databuf_create},
    obuf_desc: {hpguppi_input_databuf_create}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&thread_desc);
}

// vi: set ts=2 sw=2 et :
