#include "hpguppi_atasnap.h"

// Reset counter(s) in datablock_stats
void reset_datablock_stats(struct datablock_stats *d)
{
  d->npacket=0;
  d->ndrop=0;
}

// (Re-)initialize some or all fields of datablock_stats bi.
// d->db is set if dbout is non-NULL.
// d->block_idx is set if block_idx >= 0.
// d->block_num is always set and the stats are always reset.
// d->pkts_per_block is set of pkt_size > 0.
void init_datablock_stats(struct datablock_stats *d,
    ATA_IBV_OUT_DATABUF_T *dbout, int block_idx, int64_t block_num,
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

/* Push all blocks down a level, losing the first one */
void block_stack_push(struct datablock_stats *d, int nblock)
{
    int i;
    for (i=1; i<nblock; i++)
        memcpy(&d[i-1], &d[i], sizeof(struct datablock_stats));
}

// Update block's header info and set filled status (i.e. hand-off to downstream)
void finalize_block(struct datablock_stats *d)
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
  else if(d->pkts_per_block < d->npacket) { // Should never happen, so warn about it
    hashpipe_warn(__FUNCTION__, "block #%d has too many packets (%u > %u pkts_per_block)", d->block_idx, d->npacket, d->pkts_per_block);
  }

  sprintf(dropstat, "%d/%lu", d->ndrop, d->pkts_per_block);
  hputi4(header, "NPKT", d->npacket);
  hputi4(header, "NDROP", d->ndrop);
  hputs(header, "DROPSTAT", dropstat);
  hpguppi_databuf_set_filled(d->dbout, d->block_idx);
}

// Advance to next block in data buffer.  This new block will contain
// absolute block block_num.
//
// NB: The caller must wait for the new data block to be free after this
// function returns!
void increment_block(struct datablock_stats *d, int64_t block_num)
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
  d->packet_idx += d->pktidx_per_block;
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
void wait_for_block_free(const struct datablock_stats * d,
    hashpipe_status_t * st, const char * status_key)
{
  int rv;
  char netstat[80] = {0};
  char netbuf_status[80];
  int netbuf_full = hpguppi_databuf_total_status(d->dbout);
  sprintf(netbuf_status, "%d/%d", netbuf_full, d->dbout->header.n_block);

  hashpipe_status_lock_safe(st);
  {
    hgets(st->buf, status_key, sizeof(netstat), netstat);
    hputs(st->buf, status_key, "waitfree");
    hputs(st->buf, "NETBUFST", netbuf_status);
  }
  hashpipe_status_unlock_safe(st);

  while ((rv=hpguppi_databuf_wait_free(d->dbout, d->block_idx))
      != HASHPIPE_OK) {
    if (rv==HASHPIPE_TIMEOUT) {
    //   netbuf_full = hpguppi_databuf_total_status(d->dbout);
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

unsigned check_pkt_observability_silent(
    const struct ata_snap_obs_info * ata_oi,
    const uint64_t pkt_idx,
    const uint64_t obs_start_pktidx,
    const uint16_t feng_id,
    const uint16_t pkt_schan
  )
{
  if(pkt_idx < obs_start_pktidx){
    return PKT_OBS_IDX;
  }
  return check_pkt_observability_sans_idx(ata_oi,
                                          feng_id,
                                          pkt_schan);
}

unsigned check_pkt_observability(
    const struct ata_snap_obs_info * ata_oi,
    const uint64_t pkt_idx,
    const uint64_t obs_start_pktidx,
    const uint16_t feng_id,
    const uint16_t pkt_schan
  )
{
  unsigned obs_code = check_pkt_observability_silent(
    ata_oi,
    pkt_idx,
    obs_start_pktidx,
    feng_id,
    pkt_schan
  );
  switch(obs_code){
    case PKT_OBS_IDX:
      // hashpipe_error(__FUNCTION__, "pkt_idx (%lu) < (%lu) obs_start_pktidx", pkt_idx, obs_start_pktidx);
      break;
    case PKT_OBS_FENG:
      hashpipe_warn(__FUNCTION__, "feng_id (%u) >= (%u) ata_oi->nants", feng_id, ata_oi->nants);
      break;
    case PKT_OBS_SCHAN:
      hashpipe_warn(__FUNCTION__, "pkt_schan (%d) < (%d) ata_oi->schan", pkt_schan, ata_oi->schan);
      break;
    case PKT_OBS_NCHAN:
      hashpipe_warn(__FUNCTION__, "pkt_chans [%d-%d] <> [%d-%d] ata_oi->chans",
              pkt_schan, pkt_schan + ata_oi->pkt_nchan, ata_oi->schan, ata_oi->schan + ata_oi->nchan);
      break;
    default:
      break;
  }
  return obs_code;
}

//  if state == RECORD && STTVALID == 0 
//    STTVALID=1
//    calculate and store STT_IMJD, STT_SMJD
//  else if STTVALID != 0
//    STTVALID=0
//  endif
uint32_t update_stt_status_keys( hashpipe_status_t *st,
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
      sttvalid = 1;

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
      sttvalid = 0;
    }
    hputu4(st->buf, "STTVALID", sttvalid);
  }
  hashpipe_status_unlock_safe(st);
  return sttvalid;
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

// Check the given block's start/stop values against the observation's OBSSTART/OBSSTOP
// values (given as PKT sequence numbers). Logic goes something like this:
//   if OBSSTART <= block_start_or_stop < OBSSTOP
//     return RECORD
//   else if block_stop <= OBSSTART
//     return ARMED
//   else
//     return IDLE
//   endif
enum run_states state_from_block_start_stop(const uint64_t obs_start_pktidx, const uint64_t obs_stop_pktidx,
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

int ata_snap_obs_info_read(hashpipe_status_t *st, struct ata_snap_obs_info *obs_info)
{
  enum obs_info_validity validity = OBS_UNKNOWN; // No persistence of validity
  ata_snap_obs_info_read_with_validity(st, obs_info, &validity);
  return validity >= OBS_SEEMS_VALID;
}

// This function overwrites the status values with those of the provided ata_snap_obs_info
// which primarily serves to revert any changes that were made mid-observation
int ata_snap_obs_info_write(hashpipe_status_t *st, struct ata_snap_obs_info *obs_info)
{
  enum obs_info_validity validity = OBS_SEEMS_VALID;
  // If obs_info is valid
  if(!ata_snap_obs_info_valid(*obs_info)) {
    validity = OBS_INVALID;
  }
  ata_snap_obs_info_write_with_validity(st, obs_info, validity);
  
  return validity >= OBS_SEEMS_VALID;
}


/**
 * Updates @param validity and returns 1 if the info changed, else 0. 
 * 
 * Commented out are fields that are static due to the SNAP's implementation.
 */
char ata_snap_obs_info_read_with_validity(hashpipe_status_t *st, struct ata_snap_obs_info *obs_info, enum obs_info_validity *validity)
{
  struct ata_snap_obs_info obs_info_old = {0};
  memcpy(&obs_info_old, obs_info, sizeof(struct ata_snap_obs_info));

  // Get any obs info from status buffer, store values
  hashpipe_status_lock_safe(st);
  {
    hgetu4(st->buf, "BLOCSIZE", &obs_info->eff_block_size);
    // Read (no change if not present)
    hgetu4(st->buf, "FENCHAN",  &obs_info->fenchan);
    hgetu4(st->buf, "NANTS",    &obs_info->nants);
    hgetu4(st->buf, "OBSNCHAN", &obs_info->nchan);
    hgetu4(st->buf, "NPOL",     &obs_info->pkt_npol);
    hgetu4(st->buf, "NBITS",    &obs_info->time_nbits);
    hgetu4(st->buf, "PKTNTIME", &obs_info->pkt_ntime);
    hgetu4(st->buf, "PKTNCHAN", &obs_info->pkt_nchan);
    hgeti4(st->buf, "SCHAN",    &obs_info->schan);
    hgetr8(st->buf, "OBSBW",    &obs_info->obs_bw);
  }
  hashpipe_status_unlock_safe(st);
  obs_info->nchan /= obs_info->nants;

  // if no change in obs_info
  if (*validity != OBS_UNKNOWN &&
      obs_info_old.fenchan == obs_info->fenchan &&
      obs_info_old.nants == obs_info->nants &&
      obs_info_old.nchan == obs_info->nchan &&
      obs_info_old.pkt_npol == obs_info->pkt_npol &&
      obs_info_old.time_nbits == obs_info->time_nbits &&
      obs_info_old.pkt_ntime == obs_info->pkt_ntime &&
      obs_info_old.pkt_nchan == obs_info->pkt_nchan &&
      obs_info_old.schan == obs_info->schan &&
      obs_info_old.obs_bw == obs_info->obs_bw)
  {
    // then obs will continue to be invalid or valid
    return 0;
  }
  else if(ata_snap_obs_info_valid(*obs_info)) { // if change in obs_info and obs_info seems valid
    ata_snap_populate_block_related_fields(obs_info->eff_block_size, obs_info);
    *validity = OBS_SEEMS_VALID;
    return 1;
  } else { // if change in obs_info and obs_info seems invalid
    *validity = OBS_INVALID;
    return 1;
  }
}

// This function overwrites the status values with those of the provided ata_snap_obs_info
// which primarily serves to revert any changes that were made mid-observation
void ata_snap_obs_info_write_with_validity(hashpipe_status_t *st, struct ata_snap_obs_info *obs_info, enum obs_info_validity obs_info_valid)
{
  uint32_t obsnchan = 0;

  // Get any obs info from status buffer, store values
  hashpipe_status_lock_safe(st);
  {
    // If obs_info is valid
    if(obs_info_valid >= OBS_SEEMS_VALID) {
      obsnchan = ata_snap_obsnchan(*obs_info);
      hputs(st->buf, "OBSINFO", "VALID");
    } else {
      obsnchan = 1;
      hputs(st->buf, "OBSINFO", "INVALID");
    }

    // Write (store default/invalid values if not present)
    hputu4(st->buf, "OBSNCHAN", obsnchan);
    hputu4(st->buf, "FENCHAN",  obs_info->fenchan);
    hputu4(st->buf, "NANTS",    obs_info->nants);
    hputu4(st->buf, "NCHAN",    obs_info->nchan);
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
}

// Align the blk0_start index to be integer (preferably positive) 
// multiples of pktidx_per_block away from obsstart.
char align_blk0_with_obsstart(uint64_t * blk0_start_pktidx, uint32_t obsstart, uint32_t pktidx_per_block){
  //
  // blk0
  //  |  blk
  //  |   |  blk
  //  |   |   | observation
  //  |   |   |  |
  //  0---|---|-->____
  //          |vv|
  //    /offset/
  //  |^^|
  //  0++0---|---|>____
  //
  //  observation
  //      | blk0 
  //      |  |  blk
  //      |  |   |  blk
  //      |  |   |   |
  //      >__0___|___|
  //      |^^|offset
  //      0___|___...

  uint32_t blk_obsstart_alignment_offset = (obsstart - *blk0_start_pktidx)%pktidx_per_block;

  if(blk_obsstart_alignment_offset != 0){
    // Subtract rather, so that the offset motion is more inclusive rather than exclusive
    // (particularly for the case of blk0 > obsstart)
    if(*blk0_start_pktidx > pktidx_per_block){
      blk_obsstart_alignment_offset = pktidx_per_block - blk_obsstart_alignment_offset;
      *blk0_start_pktidx -= blk_obsstart_alignment_offset;
    }
    else{
      *blk0_start_pktidx += blk_obsstart_alignment_offset;
    }
    hashpipe_info(__FUNCTION__,
        "working blocks reinit to align pktstart to obsstart\n\t\toffset of -%lu",
        blk_obsstart_alignment_offset);
    return 1;
  }
  return 0;
}