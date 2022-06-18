// guppiraw_ibv_packetize_thread
//
// Ingest the GUPPI-RAW files with stem as per `RAWSTEM` key-value,
// generate packet-payloads emulating that of the ibverbs_packet_thread,
// and pass that buffer along. 
//
// Particularly used for offline testing, where a GUPPI-RAW file provides
// packet data.

#define _GNU_SOURCE 1
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <sys/types.h>
#include <fcntl.h>

#include "hashpipe.h"
#include "hpguppi_databuf.h"
#include "hpguppi_time.h"
#include "hpguppi_atasnap.h"
#include "hpguppi_util.h"
#include "hpguppi_pktbuf.h"
#include "hpguppi_rawspec.h"

// Parses the ibvpktsz string for chunk sizes and initializes db's pktbuf_info
// accordingly.  Returns 0 on success or -1 on error.
static
int
parse_ibvpktsz(struct hpguppi_pktbuf_info *pktbuf_info, char * ibvpktsz, size_t blocksize)
{
  int i;
  char * p;
  uint32_t nchunks = 0;
  size_t pkt_size = 0;
  size_t slot_size = 0;

  if(!ibvpktsz) {
    return -1;
  }

  // Look for commas
  while(nchunks < MAX_CHUNKS && (p = strchr(ibvpktsz, ','))) {
    // Replace comma with nul
    *p = '\0';
    // Parse chuck size
    pktbuf_info->chunks[nchunks].chunk_size = strtoul(ibvpktsz, NULL, 0);
    // Replace nul with comma
    *p = ',';
    // If chunk_size is 0, return error
    if(pktbuf_info->chunks[nchunks].chunk_size == 0) {
      hashpipe_error("IBVPKTSZ", "chunk size must be non-zero");
      return -1;
    }
    // Increment nchunks
    nchunks++;
    // Advance ibvpktsz to character beyond p
    ibvpktsz = p+1;
  }

  // If nchunks is less than MAX_CHUNKS and ibvpktsz[0] is not nul
  if(nchunks < MAX_CHUNKS && *ibvpktsz) {
    // If more commas remain, too many chunks!
    if(strchr(ibvpktsz, ',')) {
      hashpipe_error("IBVPKTSZ", "too many chunks");
      return -1;
    }
    // Parse final chunk size
    pktbuf_info->chunks[nchunks].chunk_size = strtoul(ibvpktsz, NULL, 0);
    // Increment nchunks
    nchunks++;
  } else if(nchunks == MAX_CHUNKS && *ibvpktsz) {
    // Too many chunks
    hashpipe_error("IBVPKTSZ", "too many chunks");
    return -1;
  }

  // Calculate remaining fields
  for(i=0; i<nchunks; i++) {
    pktbuf_info->chunks[i].chunk_aligned_size = pktbuf_info->chunks[i].chunk_size +
      ((-pktbuf_info->chunks[i].chunk_size) % PKT_ALIGNMENT_SIZE);
    pktbuf_info->chunks[i].chunk_offset = slot_size;
    // Accumulate pkt_size and slot_size
    pkt_size += pktbuf_info->chunks[i].chunk_size;
    slot_size += pktbuf_info->chunks[i].chunk_aligned_size;
  }

  // Store final values
  pktbuf_info->num_chunks = nchunks;
  pktbuf_info->pkt_size = pkt_size;
  pktbuf_info->slot_size = slot_size;
  pktbuf_info->slots_per_block = blocksize / slot_size;

  pktbuf_info->slots_per_block = (pktbuf_info->slots_per_block/8) * 8;

  return 0;
}

void _null_terminate_key(char *key){
  while(*key != '=' && *(key++) != ' ');
  *key = '\0';  
}

void _null_terminate_last_quote(char *str){
  str += 72;
  while(*(str--) != '\'');
  *str = '\0';  
}

void hput_buf(char *buf, char *guppiheader_buf) {
  char value[71] = {'\0'};

  // Loop over the 80-byte records until the "END " record
  while(strncmp(guppiheader_buf, "END ", 4) != 0) {
    value[0] = '\0';
    _null_terminate_key(guppiheader_buf);
    // check if target buf has the same key
    hgets(buf, guppiheader_buf, 71, value);

    if(value[0] != '\0') {
      hashpipe_info(
        __FUNCTION__,
        "Honouring existing `%s` value: `%s`",// (over '%s)",
        guppiheader_buf,
        value//,
        // guppiheader_buf+10
      );
    }
    else {
      if(*(guppiheader_buf+10) == '\''){ // strings start with " '"
        _null_terminate_last_quote(guppiheader_buf+10);
        hputs(buf, guppiheader_buf, guppiheader_buf+11);
      }
      else { // not a string, probably some numeric value
        strncpy(value, guppiheader_buf+9, sizeof(value));
        hputs(buf, guppiheader_buf, guppiheader_buf+9);
        // fprintf(stderr, "\tkey: `%s`, value: `%s`\n", guppiheader_buf, value);
      }
    }
    guppiheader_buf += 80;
  }
}

static int init(hashpipe_thread_args_t *args)
{
  // Local aliases to shorten access to args fields
  // Our output buffer happens to be a hpguppi_input_databuf
  hpguppi_input_databuf_t *db = (hpguppi_input_databuf_t *)args->obuf;
  hashpipe_status_t * st = &args->st;
  const char * status_key = args->thread_desc->skey;
  const char * thread_name = args->thread_desc->name;

  char guppifile_pathstem[73] = {'\0'};
  char guppifile_path[sizeof(guppifile_pathstem)+10] = {'\0'};

  uint32_t blockn_ant, blockn_freq, blockn_pol;
  uint32_t schan=0;
  uint32_t blocsize=0;
  uint32_t nbits=0;
  uint32_t pktnchan=0, pktntime=0;
  char ibvpktsz[80];
  strcpy(ibvpktsz, "42,16,8192");

  hashpipe_status_lock_safe(st);
    hgets(st->buf, "RAWSTEM", sizeof(guppifile_pathstem), guppifile_pathstem);
  hashpipe_status_unlock_safe(st);
  sprintf(guppifile_path, "%s.0000.raw", guppifile_pathstem);
  int guppifile_fd = open(guppifile_path, O_RDONLY);

  if (guppifile_fd == -1) {
    hashpipe_error(thread_name, "Failed to open `RAWSTEM` specified %s.", guppifile_path);
    return HASHPIPE_ERR_PARAM;
  }

  hashpipe_info(thread_name, "Opened %s.", guppifile_path);
  // read critical keys from GUPPIRAW header and push to status buffer
  char guppifile_header[MAX_RAW_HDR_SIZE];
  // Read header (plus some data, probably)
  read(guppifile_fd, guppifile_header, MAX_RAW_HDR_SIZE);
  close(guppifile_fd);

  hgetu4(guppifile_header, "BLOCSIZE", &blocsize);
  hgetu4(guppifile_header, "PKTNCHAN", &pktnchan);
  hgetu4(guppifile_header, "PKTNTIME", &pktntime);
  blockn_ant = 1; // default
  hgetu4(guppifile_header, "NANTS", &blockn_ant);
  hgetu4(guppifile_header, "OBSNCHAN", &blockn_freq);
  hgetu4(guppifile_header, "NPOL", &blockn_pol);
  hgetu4(guppifile_header, "NBITS", &nbits);
  hgetu4(guppifile_header, "SCHAN", &schan);
  fprintf(stderr, "Block parameters:\n");
  fprintf(stderr, "\tBLOCSIZE = %d\n", blocsize);
  fprintf(stderr, "\tPKTNCHAN = %d\n", pktnchan);
  fprintf(stderr, "\tPKTNTIME = %d\n", pktntime);
  fprintf(stderr, "\tNANTS    = %d\n", blockn_ant);
  fprintf(stderr, "\tOBSNCHAN = %d\n", blockn_freq);
  fprintf(stderr, "\tNPOL     = %d\n", blockn_pol);
  fprintf(stderr, "\tNBITS    = %d\n", nbits);
  fprintf(stderr, "\tSCHAN    = %d\n", schan);

  hashpipe_status_lock_safe(st);
  {
    // get keys that can override those values of the GUPPIRAW header
    hgetu4(st->buf, "PKTNCHAN", &pktnchan);
    hgetu4(st->buf, "PKTNTIME", &pktntime);

    // push keys
    hputu4(st->buf, "BLOCSIZE", blocsize);
    hputu4(st->buf, "PKTNCHAN", pktnchan);
    hputu4(st->buf, "PKTNTIME", pktntime);
    hputu4(st->buf, "NANTS", blockn_ant);
    hputu4(st->buf, "OBSNCHAN", blockn_freq);
    hputu4(st->buf, "NPOL", blockn_pol);
    hputu4(st->buf, "NBITS", nbits);
    hputu4(st->buf, "SCHAN", schan);

    hgets(st->buf, "IBVPKTSZ", sizeof(ibvpktsz), ibvpktsz);
    hputs(st->buf, "IBVPKTSZ", ibvpktsz);

    // Set status_key to init
    hputs(st->buf, status_key, "init");
  }
  hashpipe_status_unlock_safe(st);

  // blockn_time = blocsize / ((blockn_freq * blockn_pol * 2 * nbits)/8);
  // if(blockn_time % pktntime != 0) {
  //   hashpipe_error(thread_name, "NTIME (%d) %% (%d) PKTNTIME != 0", blockn_time, pktntime);
  //   return HASHPIPE_ERR_PARAM;
  // }

  // Get pointer to hpguppi_pktbuf_info
  struct hpguppi_pktbuf_info * pktbuf_info = hpguppi_pktbuf_info_ptr(db);
  // Parse ibvpktsz
  if(parse_ibvpktsz(pktbuf_info, ibvpktsz, blocsize)) {
    return HASHPIPE_ERR_PARAM;
  }
  
  if(pktbuf_info->chunks[2].chunk_size < (pktnchan*pktntime*blockn_pol*2*nbits)/8) {
    hashpipe_error(
      thread_name,
      "pktbuf_info->chunks[2].chunk_size (%d) < (%d*%d*%d*2*%d/8) pktnchan*pktntime*npol*nbits/8",
      pktbuf_info->chunks[2].chunk_size, pktnchan, pktntime, blockn_pol, nbits
    );
    return HASHPIPE_ERR_PARAM;
  }

  // Success!
  return HASHPIPE_OK;
}

static void * run(hashpipe_thread_args_t * args)
{
  // Local aliases to shorten access to args fields
  // Our input and output buffers happen to be a hpguppi_input_databuf
  struct hpguppi_input_databuf *dbout = (struct hpguppi_input_databuf *)args->obuf;
  hashpipe_status_t *st = &args->st;
  const char * thread_name = args->thread_desc->name;
  // const char * status_key = args->thread_desc->skey;

  /* administrative variables */  
  // Get pointer to hpguppi_pktbuf_info, setup in init
  size_t blockpkt_slot = 0;
  struct hpguppi_pktbuf_info * pktbuf_info = hpguppi_pktbuf_info_ptr(dbout);
  const struct hpguppi_pktbuf_chunk * chunks = pktbuf_info->chunks;
  const size_t slots_per_block = pktbuf_info->slots_per_block;
  uint32_t pktnchan=-1, pktntime=-1;
  
  uint32_t blockn_ant, blockn_freq, blockn_time, blockn_pol;
  uint16_t blockant_i=0, blockfreq_i=0, blocktime_i=0, blockpol_i=0;
  uint64_t pktidx=0;
  uint32_t schan=0;
  uint32_t blocsize=0;
  uint32_t nbits=0;

  struct ata_snap_payload_header pkt_header = {0};
  pkt_header.version = 42;
  pkt_header.type = 42;

  uint8_t* base_addr;
  char curblk = 0;//, rv;

  /* guppifile variables */
  char guppifile_pathstem[73] = {'\0'};
  char guppifile_path[sizeof(guppifile_pathstem)+10] = {'\0'};
  int guppifile_i = 0;
  off_t guppifile_pos;
  int guppifile_fd;
  rawspec_raw_hdr_t raw_hdr;

  /* start up */
  hashpipe_status_lock_safe(st);
    hgets(st->buf, "RAWSTEM", sizeof(guppifile_pathstem), guppifile_pathstem);
  
    hgetu4(st->buf, "BLOCSIZE", &blocsize);
    hgetu4(st->buf, "PKTNCHAN", &pktnchan);
    hgetu4(st->buf, "PKTNTIME", &pktntime);
    blockn_ant = 1; // default
    hgetu4(st->buf, "NANTS", &blockn_ant);
    hgetu4(st->buf, "OBSNCHAN", &blockn_freq);
    hgetu4(st->buf, "NPOL", &blockn_pol);
    hgetu4(st->buf, "NBITS", &nbits);
    hgetu4(st->buf, "SCHAN", &schan);
  hashpipe_status_unlock_safe(st);
  blockn_freq /= blockn_ant;
  
  if(chunks[2].chunk_size != (pktnchan*pktntime*blockn_pol*2*nbits)/8) {
    hashpipe_warn(
      thread_name,
      "Excessive pktpayload_bytesize is suboptimal: (%d) != (%d = %d*%d*%d*2*%d/8) pktnchan*pktntime*npol*nbits/8",
        chunks[2].chunk_size,
        (pktnchan*pktntime*blockn_pol*2*nbits)/8,
        pktnchan,
        pktntime,
        blockn_pol,
        nbits
    );
  }

  blockn_time = blocsize / ((blockn_ant * blockn_freq * blockn_pol * 2 * nbits)/8);
  hashpipe_info(thread_name, "NTIME: %d", blockn_time);

  if(blockn_time % pktntime != 0) {
    hashpipe_error(thread_name, "NTIME (%d) %% (%d) PKTNTIME != 0", blockn_time, pktntime);
    pthread_exit(NULL);
    return NULL;
  }

  if(blockn_freq % pktnchan != 0) {
    hashpipe_error(thread_name, "NCHAN (%d) %% (%d) PKTNCHAN != 0", blockn_freq, pktntime);
    pthread_exit(NULL);
    return NULL;
  }

  // GUPPI RAW data-block is (slowest)[NANT, FREQ, TIME, POL, complex-sample](fastest)
  // RTR
  const size_t blocktime_stride = (blockn_pol * 2 * nbits)/8;
  const size_t atomic_slice = pktntime * blocktime_stride;
  const size_t blockfreq_stride = blockn_time * blocktime_stride;
  const size_t blockant_stride = blockn_freq * blockfreq_stride;

  while(hpguppi_databuf_wait_free(dbout, curblk) == HASHPIPE_TIMEOUT && run_threads());

  sprintf(guppifile_path, "%s.%04d.raw", guppifile_pathstem, guppifile_i%10000);
  guppifile_fd = open(guppifile_path, O_RDONLY);

  while(guppifile_fd != -1) {
    hashpipe_info(thread_name, "Opened %s.", guppifile_path);

    if(guppifile_i == 0) {
      // read critical keys from GUPPIRAW header and push to status buffer
      char guppifile_header[MAX_RAW_HDR_SIZE];
      // Read header (plus some data, probably)
      read(guppifile_fd, guppifile_header, MAX_RAW_HDR_SIZE);
      lseek(guppifile_fd, 0, SEEK_SET);

      uint64_t pktstart=0, pktstop=0;
      hgetu8(guppifile_header, "PKTIDX", &pktstart);
      hgetu8(guppifile_header, "PKTSTART", &pktstart);
      hgetu8(guppifile_header, "PKTSTOP", &pktstop);
      hashpipe_info(thread_name, "PKTSTOP: %llu", pktstop);

      hashpipe_status_lock_safe(st);
        hputu8(st->buf, "PKTSTART", pktstart);
        hputu8(st->buf, "PKTSTOP", pktstop);
        hputs(st->buf,  "IBVSTAT", "running"); // spoof
      hashpipe_status_unlock_safe(st);
    }

    while(rawspec_raw_read_header(guppifile_fd, &raw_hdr) > 0) {
      // fprintf(stderr, "Block parameters:\n");
      // fprintf(stderr, "\tBLOCSIZE = %lu\n", raw_hdr.blocsize);
      // fprintf(stderr, "\tOBSNCHAN = %d\n",  raw_hdr.obsnchan);
      // fprintf(stderr, "\tNANTS    = %d\n",  raw_hdr.nants);
      // fprintf(stderr, "\tNBITS    = %d\n",  raw_hdr.nbits);
      // fprintf(stderr, "\tFLOATDATA= %d\n",  raw_hdr.float_data);
      // fprintf(stderr, "\tNPOL     = %d\n",  raw_hdr.npol);
      // fprintf(stderr, "\tOBSFREQ  = %g\n",  raw_hdr.obsfreq);
      // fprintf(stderr, "\tOBSBW    = %g\n",  raw_hdr.obsbw);
      // fprintf(stderr, "\tTBIN     = %g\n",  raw_hdr.tbin);
      
      if(blocsize != raw_hdr.blocsize) {
        hashpipe_error(thread_name, "BLOCSIZE changed during observation from %llu to %llu.", blocsize, raw_hdr.blocsize);
        guppifile_i = -2; // break file progression
        break;
      }
      if(pktidx != raw_hdr.pktidx) {
        hashpipe_error(thread_name, "PKTIDX %llu is not the expected %llu. Ignoring this.", raw_hdr.pktidx, pktidx);
      }
      
      // store block-data start
      guppifile_pos = lseek(guppifile_fd, 0, SEEK_CUR);
      do {
        base_addr = hpguppi_pktbuf_block_slot_ptr(dbout, curblk, blockpkt_slot);
        pkt_header.timestamp = __bswap_64(pktidx);
        pkt_header.chan = __bswap_16(schan + blockfreq_i);
        pkt_header.feng_id = __bswap_16(blockant_i);
        
        // assume first chunk is ethernet head
        // memcpy(base_addr + chunks[0].chunk_offset, &, chunks[0].chunk_size);
        // assume second chunk is packet header
        memcpy(base_addr + chunks[1].chunk_offset, &pkt_header, chunks[1].chunk_size);
        // assume third chunk is packet payload
        // have to stride frequency dimension to write PKTNCHAN into payload
        for(int chan_off = 0; chan_off <= pktnchan; chan_off++){
          lseek(
            guppifile_fd,
            guppifile_pos + (
              blocktime_i*blocktime_stride +
              (blockfreq_i + chan_off)*blockfreq_stride +
              blockant_i*blockant_stride
            ),
            SEEK_SET
          );
          read(
            guppifile_fd,
            base_addr + chunks[2].chunk_offset + chan_off*atomic_slice,
            atomic_slice
          );
        }

        // increment indices
        blockpkt_slot = (blockpkt_slot + 1)%slots_per_block;
        if(blockpkt_slot == 0){
          //! TODO update status while waiting

          // guppifile_pos = lseek(guppifile_fd, 0, SEEK_CUR);
          // lseek(guppifile_fd, raw_hdr.hdr_pos, SEEK_SET);
          // read(guppifile_fd, hpguppi_databuf_header(dbout, curblk), raw_hdr.hdr_size);
          // lseek(guppifile_fd, guppifile_pos, SEEK_SET);

          hashpipe_status_lock_safe(st);
            memcpy((char*)hpguppi_databuf_header(dbout, curblk), st->buf, HASHPIPE_STATUS_TOTAL_SIZE);
          hashpipe_status_unlock_safe(st);

          hpguppi_databuf_set_filled(dbout, curblk);
          curblk = (curblk + 1) %dbout->header.n_block;
          // progress buffers
          while(hpguppi_databuf_wait_free(dbout, curblk) == HASHPIPE_TIMEOUT && run_threads());
        }

        // Iterate through dimensions
        blockant_i = (blockant_i + 1) % blockn_ant;
        if (blockant_i == 0) {
          blockfreq_i = (blockfreq_i + pktnchan) % blockn_freq;
          if (blockfreq_i == 0) {
            blocktime_i = (blocktime_i + pktntime) % blockn_time;
            blockpol_i = (blockpol_i + blockn_pol) % blockn_pol; // practically no-op
            pktidx += pktntime;
          }
        }
      } while ((blockant_i | blockfreq_i | blocktime_i | blockpol_i) != 0);
      if(raw_hdr.directio) {
        blocsize = ((blocsize+511)/512) * 512;
      }
      lseek(guppifile_fd, guppifile_pos+blocsize, SEEK_SET);
    } // while read_header (for each block)

    close(guppifile_fd);
    if(++guppifile_i == 10000){
      break;
    }
    sprintf(guppifile_path, "%s.%04d.raw", guppifile_pathstem, guppifile_i%10000);
    guppifile_fd = open(guppifile_path, O_RDONLY);
  } // while guppifile_fd
  hashpipe_info(thread_name, "Could not open %s.", guppifile_path);

  hashpipe_info(thread_name, "exiting!");
  pthread_exit(NULL);

  return NULL;
}

static hashpipe_thread_desc_t thread_desc = {
    name: "guppiraw_ibv_packetize_thread",
    skey: "SYNSTAT",
    init: init,
    run:  run,
    ibuf_desc: {NULL},
    obuf_desc: {hpguppi_input_databuf_create}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(&thread_desc);
}

// vi: set ts=2 sw=2 et :
