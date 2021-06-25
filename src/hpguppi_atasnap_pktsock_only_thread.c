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

#define ELAPSED_S(start,stop) \
  ((int64_t)stop.tv_sec-start.tv_sec)

#define ELAPSED_NS(start,stop) \
  (ELAPSED_S(start,stop)*1000*1000*1000+(stop.tv_nsec-start.tv_nsec))

#define HPUT_DAQ_STATE(st, state)\
  hputs(st->buf, "DAQSTATE", state == IDLE  ? "idling" :\
                             state == ARMED ? "armed"  :\
                             "recording")

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
  enum obs_info_validity obs_info_validity = OBS_UNKNOWN;
  unsigned long pkt_blk_num, last_pkt_blk_num = 0;


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
  uint64_t npacket_total=0, ndrop_total=0, nbogus_total=0, npacket=0;

  uint64_t pkt_idx;
  uint64_t prev_obs_start_pktidx, prev_obs_stop_pktidx;
  uint64_t obs_start_pktidx = 0, obs_stop_pktidx = 0, blk0_start_pktidx=0;
  uint64_t pkt_payload_size;
  uint32_t fid_stride, time_stride;
  int32_t stream;
  uint16_t feng_id, pkt_schan;

  char waiting=-1, observing=0;
  char flag_reinit_blks=0;

  /* Time parameters */
  struct timespec ts_checked_obs_info = {0}, ts_tried_obs_info = {0}, ts_now = {0};
  const uint64_t obs_info_refresh_period_ns = 200*1000*1000;
  const uint64_t obs_info_retry_period_s = 5;
  uint64_t obs_info_refresh_elapsed_ns;
  
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
      // Catch any changes in OBSSTART/STOP
      clock_gettime(CLOCK_MONOTONIC, &ts_now);
      obs_info_refresh_elapsed_ns = ELAPSED_NS(ts_checked_obs_info, ts_now);
      if(obs_info_refresh_elapsed_ns > obs_info_refresh_period_ns){
        memcpy(&ts_checked_obs_info, &ts_now, sizeof(struct timespec));

        // write obs_info to overwrite any changes
        observing = pkt_idx >= obs_start_pktidx && pkt_idx < obs_stop_pktidx;
        if (obs_info_validity == OBS_VALID && // if obs_info is valid
            observing ){ //and observing
            ata_snap_obs_info_write_with_validity(st, &obs_info, obs_info_validity);
        }
        else {//otherwise read obs_info
          if(ata_snap_obs_info_read_with_validity(st, &obs_info, &obs_info_validity)){
            // this code executes if the obs_info has CHANGED
            // (ie at least once before valid observation)
            ata_snap_obs_info_write_with_validity(st, &obs_info, obs_info_validity);
            
            pkt_payload_size = ata_snap_pkt_payload_bytes(obs_info);
            fid_stride = obs_info.nstrm*pkt_payload_size;
            time_stride = obs_info.nants*fid_stride;

            memcpy(&ts_tried_obs_info, &ts_now, sizeof(struct timespec));
          }
          else if (ELAPSED_S(ts_tried_obs_info, ts_now) > obs_info_retry_period_s){
            memcpy(&ts_tried_obs_info, &ts_now, sizeof(struct timespec));
            obs_info_validity = OBS_SEEMS_VALID;
          }
        }

        // kill any observations if OBS_info_INVALID
        if (obs_info_validity < OBS_SEEMS_VALID && obs_stop_pktidx > 0){
          obs_start_pktidx = 0;
          obs_stop_pktidx = 0;
          hashpipe_status_lock_safe(st);
          {
            hputu8(st->buf, "PKTSTART", 0);
            hputu8(st->buf, "PKTSTOP", 0);
          }
          hashpipe_status_unlock_safe(st);
        }
        else{
          prev_obs_start_pktidx = obs_start_pktidx;
          prev_obs_stop_pktidx = obs_stop_pktidx;
          hashpipe_status_lock_safe(st);
          {
            hgetu8(st->buf, "PKTSTART", &obs_start_pktidx);
            hgetu8(st->buf, "PKTSTOP", &obs_stop_pktidx);
          }
          hashpipe_status_unlock_safe(st);
          
          if(obs_start_pktidx != prev_obs_start_pktidx){
            hashpipe_info(thread_name, "obs_start_pktidx changed %lu -> %lu", prev_obs_start_pktidx, obs_start_pktidx);
            if(observing){
              hashpipe_info(thread_name, "obs_start_pktidx change ignored while in observation.");
              obs_start_pktidx = prev_obs_start_pktidx;
            }
          }
          if(obs_stop_pktidx != prev_obs_stop_pktidx){
            hashpipe_info(thread_name, "obs_stop_pktidx changed %lu -> %lu", prev_obs_stop_pktidx, obs_stop_pktidx);
          }
          if(npacket > 0){// let the first reinitialisation of the blocks be due to packet discontinuity
            flag_reinit_blks = align_blk0_with_obsstart(&blk0_start_pktidx, obs_start_pktidx, obs_info.pktidx_per_block);
          }
        }

        npacket_total += npacket;

        hashpipe_status_lock_safe(st);
        {
          hputi8(st->buf, "NPKTS", npacket_total);
          hputi8(st->buf, "NDROP", ndrop_total);
          hputr4(st->buf, "PHYSPKPS", npacket*(1e9/obs_info_refresh_elapsed_ns));
          hputr4(st->buf, "PHYSGBPS", (npacket*obs_info.pkt_data_size)/((float) obs_info_refresh_elapsed_ns));

          // Overwrite any changes to reflect the static nature of these keys
          hputs(st->buf, "BINDHOST", p_ps_params->ifname);
          hputi4(st->buf, "BINDPORT", p_ps_params->port);
        }
        hashpipe_status_unlock_safe(st);
        npacket = 0;
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

    ata_snap_pkt = (struct ata_snap_pkt*) p_frame;
    // Get packet's sequence number
    pkt_idx =  ATA_SNAP_PKT_NUMBER(ata_snap_pkt);

    if (obs_info_validity < OBS_SEEMS_VALID){
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

    // Get packet's block number relative to the first block's starting index.
    pkt_blk_num = (pkt_idx - blk0_start_pktidx) / obs_info.pktidx_per_block;

    // Tally npacket between 1 Hz updates
    npacket += 1;
    // count ndrop only when a block is finalised, lest it is filled out of order.
    
    feng_id = ATA_SNAP_PKT_FENG_ID(ata_snap_pkt);
    pkt_schan = ATA_SNAP_PKT_CHAN(ata_snap_pkt);
    stream = (pkt_schan - obs_info.schan) / obs_info.pkt_nchan;

    // Only copy packet data and count packet if its wblk_idx is valid
    switch(check_pkt_observability_sans_idx(&obs_info, feng_id, stream, pkt_schan)){
      case PKT_OBS_OK:

        if(!flag_reinit_blks) {
          // Manage blocks based on pkt_blk_num
          if(pkt_blk_num == wblk[n_wblock-1].block_num + 1) {
            // Time to advance the blocks!!!

            datablock_header = datablock_stats_header(&wblk[0]);
            hputu8(datablock_header, "PKTIDX", wblk[0].packet_idx);
            hputu8(datablock_header, "BLKSTART", wblk[0].packet_idx);
            hputu8(datablock_header, "BLKSTOP", wblk[1].packet_idx);
            // Finalize first working block
            finalize_block(wblk);
            // Update ndrop counter
            ndrop_total += wblk->ndrop;
            // hashpipe_info(thread_name, "Block dropped %d packets.", wblk->ndrop);
            // Shift working blocks
            block_stack_push(wblk, n_wblock);
            // Increment last working block
            increment_block(&wblk[n_wblock-1], pkt_blk_num);
            // Wait for new databuf data block to be free
            wait_for_block_free(&wblk[n_wblock-1], st, status_key);
          }
          else if (pkt_idx < blk0_start_pktidx && blk0_start_pktidx - pkt_idx < obs_info.pktidx_per_block){
            hashpipe_info(thread_name, "pkt_idx (%lu) < (%lu) blk0_start_pktidx", pkt_idx, blk0_start_pktidx);
            hashpipe_pktsock_release_frame(p_frame);
            continue;
          }
          else if(pkt_blk_num + 1 < wblk[0].block_num //TODO dont use pkt_blk_num due to underflow
              || pkt_blk_num > wblk[n_wblock-1].block_num + 1
            ) {
            flag_reinit_blks = 1;
            blk0_start_pktidx = pkt_idx;
            align_blk0_with_obsstart(&blk0_start_pktidx, obs_start_pktidx, obs_info.pktidx_per_block);
            // Should only happen when seeing first packet when obs_info is valid
            // warn in case it happens in other scenarios
            hashpipe_warn(thread_name,
                "working blocks reinit due to packet index out of working range\n\t\t(PKTIDX %lu) [%ld, %ld  <> %lu]",
                pkt_idx, wblk[0].block_num - 1, wblk[n_wblock-1].block_num + 1, pkt_blk_num);
          }
        }

        if(flag_reinit_blks) { // Reinitialise working blocks
          flag_reinit_blks = 0;
          // Re-init working blocks for block number of current packet's block,
          // and clear their data buffers
          pkt_blk_num = (pkt_idx - blk0_start_pktidx) / obs_info.pktidx_per_block;

          for(wblk_idx=0; wblk_idx<n_wblock; wblk_idx++) {
            wblk[wblk_idx].pktidx_per_block = obs_info.pktidx_per_block;
            init_datablock_stats(wblk+wblk_idx, NULL, -1,
                pkt_blk_num+wblk_idx,
                obs_info.pkt_per_block);
            wblk[wblk_idx].packet_idx = blk0_start_pktidx + wblk[wblk_idx].block_num * obs_info.pktidx_per_block;

            // also update the working blocks' headers
            datablock_header = datablock_stats_header(wblk+wblk_idx);
            hashpipe_status_lock_safe(st);
              memcpy(datablock_header, st->buf, HASHPIPE_STATUS_TOTAL_SIZE);
            hashpipe_status_unlock_safe(st);
          }

          // trigger noteable difference in last_pkt_blk_num 
          last_pkt_blk_num = pkt_blk_num + n_wblock + 1;
        }
        
        // Update PKTIDX in status buffer if it is a new pkt_blk_num
        if(pkt_blk_num > last_pkt_blk_num || pkt_blk_num + n_wblock < last_pkt_blk_num){
          last_pkt_blk_num = pkt_blk_num;

          hashpipe_status_lock_safe(st);
            hputu8(st->buf, "BLKIDX", pkt_idx/obs_info.pktidx_per_block);
            hputu8(st->buf, "PKTIDX", pkt_idx);
            hputu8(st->buf, "BLKSTART", pkt_idx);
            hputu8(st->buf, "BLKSTOP", pkt_idx + obs_info.pktidx_per_block);
          hashpipe_status_unlock_safe(st);
        }
    
        // Once we get here, compute the index of the working block corresponding
        // to this packet.  The computed index may not correspond to a valid
        // working block!
        wblk_idx = pkt_blk_num - wblk[0].block_num;

        if(0 <= wblk_idx && wblk_idx < n_wblock) {
          // Copy packet data to data buffer of working block
          COPY_PACKET_DATA_TO_DATABUF(((struct datablock_stats*) wblk+wblk_idx),
              ata_snap_pkt->payload, (pkt_idx - blk0_start_pktidx)%obs_info.pktidx_per_block,
              feng_id, stream, pkt_schan,
              fid_stride, time_stride, pkt_payload_size, 16);//obs_info.pkt_ntime);
          // Count packet for block and for processing stats
          wblk[wblk_idx].npacket++;
        }
        else{
          // Only happens if a packet arrives late, consider n_wblock++
          hashpipe_error(thread_name, "Packet ignored: determined wblk_idx = %d", wblk_idx);
        }

        obs_info_validity = OBS_VALID;
        break;
      case PKT_OBS_FENG:
        obs_info_validity = OBS_INVALID_FENG;
        hashpipe_error(thread_name, 
          "Packet ignored: PKT_OBS_FENG\n\tfeng_id (%u) >= (%u) obs_info.nants",
          feng_id, obs_info.nants
        );
        hashpipe_status_lock_safe(st);
          hputs(st->buf, "OBSINFO", "INVALID FENG");
        hashpipe_status_unlock_safe(st);
        break;
      case PKT_OBS_SCHAN:
        obs_info_validity = OBS_INVALID_SCHAN;
        hashpipe_error(thread_name, 
          "Packet ignored: PKT_OBS_SCHAN\n\tpkt_schan (%d) < (%d) obs_info.schan",
          pkt_schan, obs_info.schan
        );
        hashpipe_status_lock_safe(st);
          hputs(st->buf, "OBSINFO", "INVALID SCHAN");
        hashpipe_status_unlock_safe(st);
        break;
      case PKT_OBS_STREAM:
        obs_info_validity = OBS_INVALID_STREAM;
        hashpipe_error(thread_name, 
          "Packet ignored: PKT_OBS_STREAM\n\tstream (%d) >= (%d) obs_info.nstrm",
          stream, obs_info.nstrm
        );
        hashpipe_status_lock_safe(st);
          hputs(st->buf, "OBSINFO", "INVALID NSTRM");
        hashpipe_status_unlock_safe(st);
        break;
      default:
        obs_info_validity = OBS_INVALID;
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
  name: "hpguppi_atasnap_pktsock_only_thread",
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


// vi: set ts=4 sw=2 et :
