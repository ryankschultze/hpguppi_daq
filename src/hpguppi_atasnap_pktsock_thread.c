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
        if(ELAPSED_NS(ts_checked_obs_startstop, ts_now) > 200*1000*1000){
          memcpy(&ts_checked_obs_startstop, &ts_now, sizeof(struct timespec));
          hashpipe_status_lock_safe(st);
          {
            hgetu8(st->buf, "PKTSTART", &obs_start_pktidx);
            hgetu8(st->buf, "PKTSTOP", &obs_stop_pktidx);
          }
          hashpipe_status_unlock_safe(st);
        }
          
        switch(state_from_start_stop(pkt_seq_num, obs_start_pktidx, obs_stop_pktidx)){
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
            hputu8(datablock_header, "BLKSTART", first_pkt_seq_num + wblk_idx * obs_info.pktidx_per_block);
            hputu8(datablock_header, "BLKSTOP", first_pkt_seq_num + (wblk_idx + 1) * obs_info.pktidx_per_block);
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
          hputu8(st->buf, "BLKSTART", pkt_seq_num);
          hputu8(st->buf, "BLKSTOP", pkt_seq_num + obs_info.pktidx_per_block);
          hashpipe_status_unlock_safe(st);
        }

        // Tally npacket_totals
        npacket_total += 1;
        // count ndrop only when a block is finalised, lest it is filled out of order.
        
        if (state != RECORD){
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
        if(pkt_blk_num == wblk[n_wblock-1].block_num + 1) {
            // Time to advance the blocks!!!
            // fprintf(stderr, "\nFinalising Block: %ld", wblk[0].block_num);
            clock_gettime(CLOCK_MONOTONIC, &ts_stop_block);

            // If the block to be finalised is out of observation range then flag_obs_end
            if(state_from_start_stop(first_pkt_seq_num + wblk[0].block_num * obs_info.pktidx_per_block, obs_start_pktidx, obs_stop_pktidx) == IDLE &&
              state_from_start_stop(first_pkt_seq_num + (wblk[0].block_num + 1) * obs_info.pktidx_per_block, obs_start_pktidx, obs_stop_pktidx) == IDLE){
                  flag_obs_end = 1;
            }

            // Finalize first working block
            datablock_header = datablock_stats_header(&wblk[0]);
            
            // update PKTIDX triggering rawdisk to close fd when this last block gets finalised.
            hputu8(datablock_header, "PKTIDX", flag_obs_end ? pkt_seq_num : first_pkt_seq_num + wblk[0].block_num * obs_info.pktidx_per_block);
            hputu8(datablock_header, "BLKSTART", first_pkt_seq_num + wblk[0].block_num * obs_info.pktidx_per_block);
            hputu8(datablock_header, "BLKSTOP", first_pkt_seq_num + (wblk[0].block_num + 1) * obs_info.pktidx_per_block);
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
              hashpipe_status_lock_safe(st);
                memcpy(datablock_header, st->buf, HASHPIPE_STATUS_TOTAL_SIZE);
              hashpipe_status_unlock_safe(st);
              hputu8(datablock_header, "PKTIDX", pkt_seq_num + wblk_idx * obs_info.pktidx_per_block);
              hputu8(datablock_header, "BLKSTART", pkt_seq_num + wblk_idx * obs_info.pktidx_per_block);
              hputu8(datablock_header, "BLKSTOP", pkt_seq_num + (wblk_idx + 1) * obs_info.pktidx_per_block);
            }
            // Immediately update STT keys to ensure that the rawdisk thread uses a new stem
            update_stt_status_keys(st, state, obs_start_pktidx, mjd);
        }

        // Check observation state
        if(flag_obs_end){// transitioning RECORD->IDLE
          flag_obs_end = 0;
          first_pkt_seq_num = 0;
          state = IDLE;
          flag_state_update = 1;
          update_stt_status_keys(st, state, pkt_seq_num, mjd);
          
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
                  ata_snap_pkt->payload, pkt_obs_relative_idx%obs_info.pktidx_per_block,
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
