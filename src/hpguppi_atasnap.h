// # hpguppi_atasnap.h - Definitions for ATA SNAP F Engine packets
//
// ## The ATA SNAP F Engine
//
// The ATA SNAP F Engine digitizes two outputs of the IF output from the ATA
// RFCBs (Radio Frequency Converter Boards) at 1800 MHz for a Nyquist bandwidth
// of 900 MHz.  The IF output of the RFCBs occupies ~600 MHz between 0 and 900
// MHz in the first Nyquist zone, so no spectral flips happen due to the
// digitization.  The two outputs typically come from two polarizations of the
// same antenna, but the gateware imposes no such hard constraints on the
// signal routing.
//
// These data are then channelized into 4096 frequency channels that are
// 900/4096 MHz (~220 kHz) wide.  The F engine can run in either "spectromweter
// mode" or "voltage mode".  In either mode, the output data is packetized and
// transmitted over a 10 GbE network interface to a user specified destination
// IP address (or addresses), which can be unitcast or multicast.  Currently,
// only voltage mode is supported by the hpguppi_atasnap_thread.
//
// In voltage mode, the complex outputs of the polyphase filter bank (PFB)
// channelizer are scaled and requantized to 4 bits for each real/imaginary
// component of the complex voltage.  With 900 MHz of bandwidth, a 4096 channel
// spectrum is produced every 4096/900e6 seconds (~4.55 usec).  At 8 bits per
// sample, the total data rate for both polarizations is 900e6*8*2 = 14.4 Gbps.
// This exceeds the data rate that can be output over a single 10 GbE link so
// only a subset of channels is sent out.  Groups of 256 contiguous channels
// are sent out over the 10 GbE network interface.  Each group is can be sent
// to an independent destination.  Each group of 256 channels constitutes 1/16
// of the data, so the data rate per group is 14.4/16 Gbps == 900 Mbps.  Up to
// 10 groups (i.e. 10/16 of the band) can be output over the 10 GbE connection
// so the maximum bandwidth available that can be output in voltage mode is 900
// MHz * 10/16 == 562.5 MHz, but for not yet understood reasons, the gateware
// only supports a maximum of 8 groups, so the maximum output bandwidth is
// limited to 450 MHz, resulting in a maximum output data rate of 7.2 Gbps (not
// counting Ethernet overhead).
//
// The channel at which the output bandwidth starts is a user selectable
// multiple of 16.
//
// Each voltage mode packet contains 2 polarizations, 256 frequency channels,
// and 16 time samples.  Since each sample is 1 byte, this leads to an 8192
// byte payload.  This payload is preceeded by an 128 bit header with the
// following format:
//
// 8  bit version
// 8  bit type
// 16 bit n_chans
// 16 bit chan
// 16 bit feng_id
// 64 bit timestamp

//
// The payload is ordered with the polarization index (0 to 1) changing
// fastest, then sample number (0 to 15), then channel number (0 to 255)
// changing the slowest.
//
// The data in a packet spans 16 spectra of 900 MHz bandwidth channelized to
// 4096 frequency channels.  This leads to a packet rate of 900e6/4096/16 ==
// ~13733 packet per second which is a packet interval of ~72.82 microseconds
// per packet.
//
// ## GUPPI RAW Blocks and Shared Memory Blocks
//
// The GUPPI RAW format stores data in fixed sized blocks.  The hpguppi network
// threads for other receivers use a shared memory block size of 128 MiB (i.e.
// 128*1024*1024).  To facilitate compatibility with these other network
// threads, the hpguppi_atasnap_thread also uses a shared memory block size
// of 128 MiB, but it creates GUPPI RAW blocks that are that size or smaller,
// depending on the nuber of antennas in use in a given observation.
//
// The GUPPI RAW block size must also be an integer multiple of the packet size
// to avoid compilcations arising from splitting a packet across block
// boundaries.
//
// Another constraint on block size is that the real-time spectroscopy code
// (i.e. rawspec) only works on an integer number of blocks.  This means that
// the highest spectral resolution desired must use an integer multiple of
// blocks.  The spectroscopy code also places constraints on the resolution of
// other (lower resolution) products computed at the same time as well as their
// integration times, but those constraints do not affect the block sizing.
//
// ## Block Size and Spectral Resolution
//
// Because the channel width of the voltage data is 900/4096 MHz, there is no
// possibility for a fine channelization that would result in an integer number
// of spectra for a small integer number of seconds.  This makes the choice of
// block size somewhat simpler in that we do not have the opportunity to choose
// between a powers of two number of fine channels per coarse channel and a
// non-power of two.  Given that the desired finest spectral resolution is ~1
// Hz, the power of two number of fine channels per coarse channel that get the
// closest to that is one of:
//
//     900/4096 MHz / 2**17 channels == 1.676 Hz/channel, 0.597 sec/spectrum
//     900/4096 MHz / 2**18 channels == 0.838 Hz/channel, 1.193 sec/spectrum
//
// Because these are both powers of two, we can set the number of time samples
// per block to a power of two less than 2**17 and an integer number of blocks
// will produce either of these two spectral resolutions.
// 
// This leaves NANTS as the parameter that defines the GUPPI RAW block size.
// Given a Hashpipe block size of 128 MiB and packets with 8192 samples in
// them, we can store up to 16384 packets per Hashpipe block.  The actual
// number of time samples we can store per Hashpipe block is:
//
//     NTIME = 128 MiB / 2 pols / NCHAN / 2**nextpow2(NANTS)
//
// Which leads to a GUPPI RAW block size (BLOCSIZE) of:
//
//     BLOCSIZE = 2 pols * NCHAN * NTIME * NANTS
//
// ## Time representation
//
// The header of the voltage mode packets contain a 38 bit counter that
// indicates the spectrum number (relative to the last synchronizing of the
// SNAP) of the first time sample in the packet payload.  At 900e6/4096 spectra
// per second, this 38 bit counter will roll over after approximately two
// weeks (~14.48 days).
//
// GUPPI RAW uses a number of independent fields to represent time in different
// ways, but the most precise one is PXTIDX (packet index).  PKTIDX is a
// monotonically increasing counter that has a direct relationship to elapsed
// time since the counter was last reset.  By knowing this relationship and the
// time at which the counter was reset, the absolute time corresponding to a
// given PKTIDX can be calculated.  In practice, conversion to absolute time is
// rarely performed.  More often, absolute times (scan start, scan stop, etc)
// are converted into PKTIDX values.  PKTIDX is used to reassemble the packets
// into a continuous sequence of data and to control the start/stop of
// recording based on absolute times that have been converted to values in the
// same timebase as PKTIDX.  PKTIDX is typically stored as a 64 bit unsigned
// integer, but it is generally safe and sometimes more convenient to store it
// as a 64 bit signed integer and treat it as a 63 bit unsigned integer.
//
// The 38 bit "sample number" in the packet headers could be used natively as
// PKTIDX, but there are a few straightforward ways in which it can be
// improved.  The first is to right shift the value by 4 bits since this
// counter is always a multiple of 16 (TODO Verify whether that's actually
// true!).  The other is to add a "rollover" count in the bits above the
// topmost bit of the sample number so that PKTIDX will not roll over when
// sample number does.  Doing these two steps will prevent PKTIDX from ever
// rolling over (practically speaking, that is; if you wait around 2 billion
// weeks or so you would get to see PKTIDX rollover).  Keeping track of the
// rollover count for the "sample number" counter requires knowledge of the
// time at which the counter was reset.  This value is generally stored
// as the number of seconds since the UIX epoch (1970-01-01T00:00:00 UTC)
// when the FPGA is synchronized to a 1 PPS pulse.  The SANP initializtion code
// actually stores this value in a register on the SANP board itself, but
// ideally this value will also be stored in the Hashpipe status buffer as
// SYNCTIME.  The status buffer gets written out as the header of each GUPPI
// RAW block so the SYNCTIME and PKTIDX values will be avilable in the output
// products.  Another useful status buffer field relating to timekeeping is
// "PIPERBLK", which is short for "PktIdx PER BLocK".  PIPERBLK is the PKTIDX
// step from one block to the next.  Its value depends on NTIME.
//
// ## IBVPKTSZ and voltage mode packet layout
//
// The ATA SNAP voltage mode packets consist of 42 bytes of Ethernet+IP+UDP
// headers, 16 bytes of application header, and 8192 bytes of payload data.
// When using the ibvpkt thread (currently known as
// "hpguppi_ibverbs_pkt_thread"), the "IBVPKYSZ" status keyword can be
// specified as "42,16,8192".

#ifndef _HPGUPPI_ATASNAP_H_
#define _HPGUPPI_ATASNAP_H_

#define XGPU_BLOCK_NANTS 32 // TFP transposition related definition

#define ATA_PAYLOAD_TRANSPOSE_FTP 0
#define ATA_PAYLOAD_TRANSPOSE_TFP 1
#define ATA_PAYLOAD_TRANSPOSE_TFP_DP4A 2

#define ATA_PAYLOAD_TRANSPOSE ATA_PAYLOAD_TRANSPOSE_FTP
// #define ATA_PACKET_PAYLOAD_DIRECT_COPY // define to use assignment copy in place of memcpy

#if ATA_PAYLOAD_TRANSPOSE == ATA_PAYLOAD_TRANSPOSE_FTP
#define COPY_PACKET_PAYLOAD_FORLOOP COPY_PACKET_DATA_TO_FTP_DATABUF_FORLOOP
#ifdef ATA_PACKET_PAYLOAD_DIRECT_COPY
  #define PKT_PAYLOAD_CP_T PKT_DCP_FTP_T
#else
  #define PKT_PAYLOAD_CP_T uint8_t
#endif

#elif ATA_PAYLOAD_TRANSPOSE == ATA_PAYLOAD_TRANSPOSE_TFP
#define COPY_PACKET_PAYLOAD_FORLOOP COPY_PACKET_DATA_TO_TFP_DATABUF_FORLOOP
#ifdef ATA_PACKET_PAYLOAD_DIRECT_COPY
  #define PKT_PAYLOAD_CP_T PKT_DCP_TFP_T
#else
  #define PKT_PAYLOAD_CP_T uint8_t
#endif

#elif ATA_PAYLOAD_TRANSPOSE == ATA_PAYLOAD_TRANSPOSE_TFP_DP4A
#define COPY_PACKET_PAYLOAD_FORLOOP COPY_PACKET_DATA_TO_TFP_DP4A_DATABUF_FORLOOP
#define PKT_PAYLOAD_CP_T PKT_DCP_TFP_DP4A_T // always a direct copy

#endif // ATA_PAYLOAD_TRANPOSE == ***

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#include "hashpipe.h"
#include "hpguppi_time.h"
#include "hpguppi_util.h"
#include "hashpipe_packet.h"

#if ATA_PAYLOAD_TRANSPOSE == ATA_PAYLOAD_TRANSPOSE_FTP
#include "hpguppi_databuf.h"
#define ATA_IBV_OUT_DATABUF_CREATE hpguppi_input_databuf_create
#define ATA_IBV_OUT_DATABUF_T struct hpguppi_input_databuf
#else
#include "hpguppi_xgpu_databuf.h"
#define ATA_IBV_OUT_DATABUF_CREATE hpguppi_input_xgpu_databuf_create
#define ATA_IBV_OUT_DATABUF_T struct hpguppi_input_xgpu_databuf
#endif

hashpipe_databuf_t *hpguppi_ata_ibv_output_databuf_create(int instance_id, int databuf_id);

static inline ATA_IBV_OUT_DATABUF_T *hpguppi_ata_ibv_output_databuf_attach(int instance_id, int databuf_id)
{
    return (ATA_IBV_OUT_DATABUF_T *)hashpipe_databuf_attach(instance_id, databuf_id);
}

static inline int hpguppi_ata_ibv_output_databuf_detach(ATA_IBV_OUT_DATABUF_T *d)
{
    return hashpipe_databuf_detach((hashpipe_databuf_t *)d);
}

static inline void hpguppi_ata_ibv_output_databuf_clear(ATA_IBV_OUT_DATABUF_T *d)
{
    hashpipe_databuf_clear((hashpipe_databuf_t *)d);
}

static inline int hpguppi_ata_ibv_output_databuf_block_status(ATA_IBV_OUT_DATABUF_T *d, int block_id)
{
    return hashpipe_databuf_block_status((hashpipe_databuf_t *)d, block_id);
}

static inline int hpguppi_ata_ibv_output_databuf_total_status(ATA_IBV_OUT_DATABUF_T *d)
{
    return hashpipe_databuf_total_status((hashpipe_databuf_t *)d);
}

static inline int hpguppi_ata_ibv_output_databuf_wait_free_timeout(
    ATA_IBV_OUT_DATABUF_T *d, int block_id, struct timespec *timeout)
{
    return hashpipe_databuf_wait_free_timeout((hashpipe_databuf_t *)d,
        block_id, timeout);
}

static inline int hpguppi_ata_ibv_output_databuf_wait_free(ATA_IBV_OUT_DATABUF_T *d, int block_id)
{
    return hashpipe_databuf_wait_free((hashpipe_databuf_t *)d, block_id);
}

static inline int hpguppi_ata_ibv_output_databuf_busywait_free(ATA_IBV_OUT_DATABUF_T *d, int block_id)
{
    return hashpipe_databuf_busywait_free((hashpipe_databuf_t *)d, block_id);
}

static inline int hpguppi_ata_ibv_output_databuf_wait_filled_timeout(
    ATA_IBV_OUT_DATABUF_T *d, int block_id, struct timespec *timeout)
{
    return hashpipe_databuf_wait_filled_timeout((hashpipe_databuf_t *)d,
        block_id, timeout);
}

static inline int hpguppi_ata_ibv_output_databuf_wait_filled(ATA_IBV_OUT_DATABUF_T *d, int block_id)
{
    return hashpipe_databuf_wait_filled((hashpipe_databuf_t *)d, block_id);
}

static inline int hpguppi_ata_ibv_output_databuf_busywait_filled(ATA_IBV_OUT_DATABUF_T *d, int block_id)
{
    return hashpipe_databuf_busywait_filled((hashpipe_databuf_t *)d, block_id);
}

static inline int hpguppi_ata_ibv_output_databuf_set_free(ATA_IBV_OUT_DATABUF_T *d, int block_id)
{
    return hashpipe_databuf_set_free((hashpipe_databuf_t *)d, block_id);
}

static inline int hpguppi_ata_ibv_output_databuf_set_filled(ATA_IBV_OUT_DATABUF_T *d, int block_id)
{
    return hashpipe_databuf_set_filled((hashpipe_databuf_t *)d, block_id);
}

static inline char *hpguppi_ata_ibv_output_databuf_header(ATA_IBV_OUT_DATABUF_T *d, int block_id) {
    if(block_id < 0 || d->header.n_block < block_id) {
        hashpipe_error(__FUNCTION__,
            "block_id %s out of range [0, %d)",
            block_id, d->header.n_block);
        return NULL;
    } else {
        return d->block[block_id].hdr;
    }
}

static inline char *hpguppi_ata_ibv_output_databuf_data(ATA_IBV_OUT_DATABUF_T *d, int block_id) {
    if(block_id < 0 || d->header.n_block < block_id) {
        hashpipe_error(__FUNCTION__,
            "block_id %s out of range [0, %d)",
            block_id, d->header.n_block);
        return NULL;
    } else {
        return d->block[block_id].data;
    }
}

// Define run states.  Currently three run states are defined: IDLE, LISTEN,
// and RECORD.
//
// In the LISTEN and RECORD states, the PKTIDX field is updated with the value
// from received packets.  Whenever the first PKTIDX of a block is received
// (i.e. whenever PKTIDX is a multiple of pktidx_per_block), the value
// for PKTSTART and DWELL are read from the status buffer.  PKTSTART is rounded
// down, if needed, to ensure that it is a multiple of pktidx_per_block,
// then PKTSTART is written back to the status buffer.  DWELL is interpreted as
// the number of seconds to record and is used to calculate PKTSTOP (which gets
// rounded down, if needed, to be a multiple of pktidx_per_block).
//
// The IDLE state is entered when there is no DESTIP defined in the status
// buffer or it is 0.0.0.0.  In the IDLE state, the DESTIP value in the status
// buffer is checked once per second.  If it is found to be something other
// than 0.0.0.0, the state transitions to the LISTEN state and the current
// blocks are reinitialized.
//
// To be operationally compatible with other hpguppi net threads, a "command
// FIFO" is created and read from in all states, but commands sent there are
// ignored.  State transitions are controlled entirely by DESTIP and
// PKTSTART/DWELL status buffer fields.
//
// In the LISTEN state, incoming packets are processed (i.e. stored in the net
// thread's output buffer) and full blocks are passed to the next thread.  When
// the processed PKTIDX is equal to PKTSTART the state transitions to RECORD
// and the following actions occur:
//
//   1. The MJD of the observation start time is calculated from PKTIDX,
//      SYNCTIME, and other parameters.
//
//   2. The packet stats counters are reset
//
//   3. The STT_IMDJ and STT_SMJD are updated in the status buffer
//
//   4. STTVALID is set to 1
//
// In the RECORD state, incoming packets are processed (i.e. stored in the net
// thread's output buffer) and full blocks are passed to the next thread (same
// as in the LISTEN state).  When the processed PKTIDX is greater than or equal
// to PKTSTOP the state transitions to LISTEN and STTVALID is set to 0.
//
// The PKTSTART/PKTSTOP tests are done every time the work blocks are advanced.
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
                    PKT_OBS_NCHAN
                  };
enum obs_info_validity { 
                OBS_UNKNOWN=-5,       //=-5
                OBS_INVALID_FENG=-4,  //=-4
                OBS_INVALID_SCHAN=-3, //=-3
                OBS_INVALID_NCHAN=-2, //=-2
                OBS_INVALID=-1,       //=-1
                OBS_SEEMS_VALID=0,    //=0   sign bit for validity
                OBS_VALID=1           //=1
              };

struct __attribute__ ((__packed__)) ata_snap_payload_header {
  uint8_t version;
  uint8_t type;
  uint16_t n_chans;
  uint16_t chan;
  uint16_t feng_id;
  uint64_t timestamp;
};

/* Structs/functions to more easily deal with multiple
 * active blocks being filled
 */
struct datablock_stats {
    ATA_IBV_OUT_DATABUF_T *dbout; // Pointer to overall shared mem databuf
    int block_idx;                    // Block index number in databuf
    int64_t block_num;                // Absolute block number
    uint64_t packet_idx;              // Index of first packet number in block
    int pktidx_per_block;            // Total number of packets to go in the block
    uint64_t pkts_per_block;
    int npacket;                      // Number of packets filled so far
    int ndrop;                     // Number of dropped packets so far
    uint64_t last_pkt;                // Last packet seq number written to block
};

// ATA SNAP packet with link layer header and internal padding to optimize
// alignment.  The alignment is acheived through judicious use of IB Verbs
// scatter/gather capabilities (specifically the scatter part).
struct __attribute__ ((__packed__)) ata_snap_ibv_pkt {
  struct ethhdr ethhdr; //0
  struct iphdr iphdr;   //14
  struct udphdr udphdr; //34
  uint8_t pad0[22];     //42
  struct ata_snap_payload_header snaphdr;      //64
  uint8_t pad1[48];     //80
	uint8_t payload[];    //128
};

struct __attribute__ ((__packed__)) ata_snap_pkt {
  struct ethhdr ethhdr; // headers aren't functional  // 14 bytes
  struct iphdr iphdr; // headers aren't functional    // 20 bytes
  struct udphdr udphdr; // headers aren't functional  // 8 bytes
  uint8_t pad0[66]; // Dont know the reason, but actual network payload is offset
  struct ata_snap_payload_header snaphdr;
  uint8_t payload[];//complex4 data[n_chans, 16, 2] // 4-bit real + 4-bit imaginary
};

#define ATA_SNAP_PKT_NUMBER(ata_snap_pkt)   (uint64_t)__bswap_64(ata_snap_pkt->snaphdr.timestamp)
#define ATA_SNAP_PKT_SCHAN(ata_snap_pkt)    __bswap_16(ata_snap_pkt->snaphdr.chan)
#define ATA_SNAP_PKT_NCHAN(ata_snap_pkt)    __bswap_16(ata_snap_pkt->snaphdr.n_chans)
#define ATA_SNAP_PKT_FENG_ID(ata_snap_pkt)  __bswap_16(ata_snap_pkt->snaphdr.feng_id)

// ATA SNAP header byte offset within (unpadded) packet
#define ATA_SNAP_PKT_OFFSET_HEADER \
  (sizeof(struct ethhdr) + \
   sizeof(struct iphdr ) + \
   sizeof(struct udphdr))

#define ATA_SNAP_PKT_SIZE_HEADER (sizeof(struct ata_snap_payload_header))

// ATA SNAP payload byte offset within (unpadded) packet
#define ATA_SNAP_PKT_OFFSET_PAYLOAD \
  (ATA_SNAP_PKT_OFFSET_HEADER + ATA_SNAP_PKT_SIZE_HEADER)

#define ATA_SNAP_PKT_SIZE_PAYLOAD (8192)

// The ATA SNAP design send out 8192 bytes of payload per packet.  MAX_PKT_SIZE
// is the next power of two that fits that plus non-payload bytes.
// TODO Why next power of two???
#define MAX_PKT_SIZE (16384)

struct ata_snap_pkt_info {
  uint64_t pktidx;
  uint16_t feng_id;
  uint16_t pkt_schan;
  uint16_t pkt_nchan;
};

// Parameters describing various data dimensions for an ATA SNAP observation.
// Only schan can meaningfully be zero, all other fields should be non-zero.
// Structure is valid if all fields are non-zero except for schan which is
// valid except for -1.
struct ata_snap_obs_info {
  // Total numner of F Engine channels
  uint32_t fenchan;
  // Total number of antennas in current subarray
  uint32_t nants;
  // Number of F Engine channels to be processed
  uint32_t nchan;
  // Number of polarisations per time sample
  uint32_t pkt_npol;
  // Number of bits per time sample component (real/imaginary)
  uint32_t time_nbits;
  // Number of time samples per packet
  uint32_t pkt_ntime;
  // Number of frequency channels per packet
  uint32_t pkt_nchan;
  // Starting F Engine channel number to be processed
  int32_t schan;
  // The number of bytes per packets
  uint32_t pkt_data_size;
  // The number of (effective) packets per block
  uint32_t pkt_per_block;
  // The number of packet timestamps per block
  uint32_t pktidx_per_block;
  // The effective size of the block
  uint32_t eff_block_size;
  // The Bandwidth of the observation in MHz
  double obs_bw;
};

#define ATASNAP_DEFAULT_FENCHAN         (  4096)
#define ATASNAP_DEFAULT_NCHAN           (   128)
#define ATASNAP_DEFAULT_PKTNPOL         (     2)
#define ATASNAP_DEFAULT_TIME_NBITS      (     4)
#define ATASNAP_DEFAULT_PKTNCHAN        (   256)
#define ATASNAP_DEFAULT_PKTNTIME        (    16)
#define ATASNAP_DEFAULT_PKT_SIZE        (  8208)
#define ATASNAP_DEFAULT_PKT_PER_BLK     ( 16834)
#define ATASNAP_DEFAULT_PKTIDX_PER_BLK  (262144)
#define OBS_INFO_INVALID_SCHAN          (    -1)

// Returns the largest power of two that it less than or equal to x.
// Returns 0 if x is 0.
static inline
uint32_t
prevpow2(uint32_t x)
{
  return x == 0 ? 0 : (0x80000000 >> __builtin_clz(x));
}

// Initialize all ata_snap_obs_info fields to invalid values
static inline
void
ata_snap_obs_info_init(struct ata_snap_obs_info * poi)
{
  memset(poi, 0, sizeof(struct ata_snap_obs_info));
  poi->fenchan = ATASNAP_DEFAULT_FENCHAN;
  poi->nchan = ATASNAP_DEFAULT_NCHAN;
  poi->pkt_npol = ATASNAP_DEFAULT_PKTNPOL;
  poi->time_nbits = ATASNAP_DEFAULT_TIME_NBITS;
  poi->pkt_ntime = ATASNAP_DEFAULT_PKTNTIME;
  poi->pkt_nchan = ATASNAP_DEFAULT_PKTNCHAN;
  poi->schan = OBS_INFO_INVALID_SCHAN;
  poi->pkt_data_size = ATASNAP_DEFAULT_PKT_SIZE;
  poi->pkt_per_block = ATASNAP_DEFAULT_PKT_PER_BLK;
  poi->pktidx_per_block = ATASNAP_DEFAULT_PKTIDX_PER_BLK;
}

static inline
int
ata_snap_obs_info_valid(const struct ata_snap_obs_info oi)
{
  return
    (oi.fenchan    != 0) &&
    (oi.nants      != 0) &&
    (oi.nchan      != 0) &&
    (oi.pkt_npol   != 0) &&
    (oi.time_nbits != 0) &&
    (oi.pkt_ntime  != 0) &&
    (oi.pkt_nchan  != 0) &&
    (oi.schan      != OBS_INFO_INVALID_SCHAN);
}

// For ATA SNAP, the OBSNCHAN parameter represents the total number of
// frequency channels processed.  It is the number of antennas times the number
// of channels per antenna.
static inline
uint32_t
calc_ata_snap_obsnchan(uint32_t nants, uint32_t nchan)
{
  return nants * nchan;
}

static inline
uint32_t
ata_snap_obsnchan(const struct ata_snap_obs_info oi)
{
  return calc_ata_snap_obsnchan(oi.nants, oi.nchan);
}

// This function doesn't assume each sample is 2 bytes ([4 bits real + 4 bits imag] *
// 2 pols).
static inline
uint32_t
calc_ata_snap_pkt_payload_bytes(uint32_t pkt_nchan,
                               uint32_t pkt_ntime, uint32_t pkt_npol,
                               uint32_t time_nbits)
{ // each time sample is real
  return (pkt_nchan * pkt_ntime * pkt_npol * 2 * time_nbits) / 8;
}

static inline
uint32_t
ata_snap_pkt_payload_bytes(const struct ata_snap_obs_info oi)
{ 
  return calc_ata_snap_pkt_payload_bytes(oi.pkt_nchan, oi.pkt_ntime, oi.pkt_npol, oi.time_nbits);
}

// Calculate the number of pktidx values per block.  Note that nchan is the
// total number of channels across all F engines.  Each PKTIDX corresponds to a
// set of packets that all share a common PKTIDX. It is possible that the
// number of PKTIDX values per block (i.e. packet sets per block) times the
// number of time samples per packet will not divide the data block size evenly
// (e.g. when NANTS is not a power of two).  This means that some trailing
// bytes in the data buffer will be unoccupied.  These unoccupied bytes should
// not be processed (e.g. copied to output files) so it is important to update
// the BLOCSIZE field in the status buffer accordingly whenever a new
// pktidx_per_block value is calculated.
static inline
uint32_t
calc_ata_snap_pkt_per_block(size_t block_size, uint32_t pkt_nchan,
                               uint32_t pkt_ntime, uint32_t pkt_npol,
                               uint32_t time_nbits)
{
  uint32_t pkt_idxpblk = (uint32_t)(block_size /
      calc_ata_snap_pkt_payload_bytes(pkt_nchan, pkt_ntime, pkt_npol, time_nbits));
  return pkt_idxpblk;
}

static inline
uint32_t
ata_snap_pkt_per_block(size_t block_size, const struct ata_snap_obs_info oi)
{
  return calc_ata_snap_pkt_per_block(
                  block_size, oi.pkt_nchan, oi.pkt_ntime,
                  oi.pkt_npol, oi.time_nbits);
}

static inline
uint32_t
calc_ata_snap_eff_pkt_per_block(size_t block_size, uint32_t pkt_nchan,
                               uint32_t pkt_ntime, uint32_t pkt_npol,
                               uint32_t time_nbits, uint32_t nants,
                               uint32_t nchan)
{
  // use integer division to round down to nearest multiple of (nants*nstrm)
  uint32_t total_strms = nants * (nchan/pkt_nchan);
  return total_strms *
        (calc_ata_snap_pkt_per_block(block_size, pkt_nchan, pkt_ntime, pkt_npol, time_nbits)
        / total_strms);
}

static inline
uint32_t
ata_snap_eff_pkt_per_block(size_t block_size, const struct ata_snap_obs_info oi)
{
  return calc_ata_snap_eff_pkt_per_block(
                  block_size, oi.pkt_nchan, oi.pkt_ntime,
                  oi.pkt_npol, oi.time_nbits,
                  oi.nants, oi.nchan);
}

static inline
uint32_t
calc_ata_snap_pkt_bytes(uint32_t pkt_nchan,
                        uint32_t pkt_ntime, uint32_t pkt_npol,
                        uint32_t time_nbits)
{
  return calc_ata_snap_pkt_payload_bytes(pkt_nchan, pkt_ntime, pkt_npol, time_nbits) + 16;
}

static inline
uint32_t
ata_snap_pkt_bytes(const struct ata_snap_obs_info oi)
{ 
  return calc_ata_snap_pkt_bytes(oi.pkt_nchan, oi.pkt_ntime, oi.pkt_npol, oi.time_nbits);
}

// For ATA SNAP, the PIPERBLK parameter (Packet Index Per Block)
// represents the total number of time samples in a block, per channel. It
// depends on the block size and NCHAN (and sample size, but that is assumed to
// be 2 bytes).  This calculation is a bit tricky because the effective block
// size can be less than the max block size when NCHAN and PKT_NTIME do not
// evenly divide the max block size.
static inline
uint32_t
calc_ata_snap_pktidx_per_block(size_t block_size, uint32_t pkt_nchan, uint32_t pkt_ntime,
                            uint32_t pkt_npol, uint32_t time_nbits,
                            uint32_t nants, uint32_t nchan)
{
  return (calc_ata_snap_eff_pkt_per_block(block_size, 
                                            pkt_nchan, pkt_ntime, pkt_npol, time_nbits,
                                            nants, nchan)
                              /(nants*(nchan/pkt_nchan))) * pkt_ntime;
}

static inline
uint32_t
ata_snap_pktidx_per_block(size_t block_size, const struct ata_snap_obs_info oi)
{
  return calc_ata_snap_pktidx_per_block(block_size, oi.pkt_nchan,
                                     oi.pkt_ntime, oi.pkt_npol, oi.time_nbits,
                                     oi.nants, oi.nchan);
}

// Calculate the effective block size for the given max block size, nchan, and
// pkt_ntime values.  The effective block size can be less than the max block size
// if nchan and/or pkt_ntime do not evenly divide the max block size.
static inline
uint32_t
calc_ata_snap_block_size(size_t block_size, uint32_t pkt_nchan, uint32_t pkt_ntime,
                         uint32_t pkt_npol, uint32_t time_nbits, uint32_t nants, uint32_t nchan)
{
  return    calc_ata_snap_pkt_payload_bytes(pkt_nchan, pkt_ntime,
                                            pkt_npol, time_nbits)
            * calc_ata_snap_eff_pkt_per_block(
                  block_size, pkt_nchan, pkt_ntime,
                  pkt_npol, time_nbits,
                  nants, nchan);
}

static inline
uint32_t
ata_snap_block_size(size_t block_size, const struct ata_snap_obs_info oi)
{
  return calc_ata_snap_block_size(block_size, oi.pkt_nchan,
                                     oi.pkt_ntime, oi.pkt_npol, oi.time_nbits,
                                     oi.nants, oi.nchan);
}

static inline
double
calc_ata_snap_tbin(uint32_t obsnchan, double obsbw_mhz)
{
  return ((double)obsnchan) / obsbw_mhz / 1e6;
}

static inline
double
ata_snap_tbin(const struct ata_snap_obs_info oi)
{
  return calc_ata_snap_tbin(ata_snap_obsnchan(oi), oi.obs_bw);
}


// The number of time samples per block (PIPERBLK) is desired to be a power
// of two, so that subsequent FFTs can operate on complete blocks with maximum
// efficiency (CUFFT preference).  For this to happen, both the number of time 
// samples per packet (PKT_NTIME) and the number of PKTIDX values per block
// (PIPERBLK) must be powers of two. PKT_NTIME is set by the upstream F engines,
// so we have no control of that, but we can and do ensure that PIPERBLK is a
// power of 2.
static inline
void
ata_snap_populate_block_related_fields(size_t block_size, struct ata_snap_obs_info *oi)
{
  oi->pkt_data_size = ata_snap_pkt_bytes(*oi);
  oi->pkt_per_block = ata_snap_eff_pkt_per_block(block_size, *oi);
  oi->pktidx_per_block = ata_snap_pktidx_per_block(block_size, *oi);//inherently effective 
  if(prevpow2(oi->pktidx_per_block) != oi->pktidx_per_block){
    hashpipe_warn(__FUNCTION__, "Recommended that the BLOCK_DATA_SIZE be changed to %lu in order to have power of 2 (%d) packets per block.", 
      (oi->pkt_data_size-16)*((prevpow2(oi->pktidx_per_block)/oi->pkt_ntime)*oi->nants*(oi->nchan/oi->pkt_nchan)), prevpow2(oi->pktidx_per_block));
  }
  oi->eff_block_size = oi->pkt_per_block*(oi->pkt_data_size-16);
}

// Parses a ATA SNAP F Engine packet that is in the format of a "struct
// ata_snap_ibv_pkt", stores metadata in fei, returns pointer to packet's
// payload.
static inline
const uint8_t *
ata_snap_parse_ibv_packet(const struct ata_snap_ibv_pkt *p,
    struct ata_snap_pkt_info * fei)
{
  fei->pktidx = ATA_SNAP_PKT_NUMBER(p);
  fei->feng_id = ATA_SNAP_PKT_FENG_ID(p);
  fei->pkt_schan = ATA_SNAP_PKT_SCHAN(p);
  fei->pkt_nchan = ATA_SNAP_PKT_NCHAN(p);

  return p->payload;
}

#if 0
static inline
uint64_t
atasnap_get_pktidx(struct vdifhdr * p)
{
  uint64_t pktidx = ((uint64_t)PKSUWL_PKTIDX_PER_SEC) * vdif_get_time(p);
  pktidx += vdif_get_data_frame_seq(p);
  return pktidx;
}

static inline
void
pksuwl_pktidx_to_timeval(uint64_t pktidx, struct timeval *tv)
{
  tv->tv_sec = pktidx / PKSUWL_PKTIDX_PER_SEC;
  tv->tv_usec = (pktidx % PKSUWL_PKTIDX_PER_SEC) * PKSUWL_NS_PER_PKT / 1000;
}
#endif

#if 0
// Not used, but here it is if needed someday...
static inline
void
pksuwl_pktidx_to_timespec(uint64_t pktidx, struct timespec *ts)
{
  ts->tv_sec = pktidx / PKSUWL_PKTIDX_PER_SEC;
  ts->tv_nsec = (pktidx % PKSUWL_PKTIDX_PER_SEC) * PKSUWL_NS_PER_PKT;
}
#endif

static inline void hput_obsdone(hashpipe_status_t * st, const int flag)
{
  hashpipe_status_lock_safe(st);
  {
    hputi4(st->buf, "OBSDONE", flag);
  }
  hashpipe_status_unlock_safe(st);
}

static inline void hget_obsdone(hashpipe_status_t * st, int *flag)
{
  hashpipe_status_lock_safe(st);
  {
    hgeti4(st->buf, "OBSDONE", flag);
  }
  hashpipe_status_unlock_safe(st);
}

// Returns pointer to datablock_stats's output data block
static inline char * datablock_stats_data(const struct datablock_stats *d)
{
  return hpguppi_databuf_data(d->dbout, d->block_idx);
}

// Returns pointer to datablock_stats's header
static inline char * datablock_stats_header(const struct datablock_stats *d)
{
  return hpguppi_databuf_header(d->dbout, d->block_idx);
}

void reset_datablock_stats(struct datablock_stats *d);
void init_datablock_stats(struct datablock_stats *d,
    ATA_IBV_OUT_DATABUF_T *dbout, int block_idx, int64_t block_num,
    uint64_t pkts_per_block);
void block_stack_push(struct datablock_stats *d, int nblock);
void finalize_block(struct datablock_stats *d);
void increment_block(struct datablock_stats *d, int64_t block_num);
void wait_for_block_free(const struct datablock_stats * d,
  hashpipe_status_t * st, const char * status_key);

static inline
unsigned check_pkt_observability_sans_idx(
    const struct ata_snap_obs_info * ata_oi,
    const uint16_t feng_id,
    const uint16_t pkt_schan
  )
{
  if(feng_id >= ata_oi->nants){
    return PKT_OBS_FENG;
  }
  if(pkt_schan < ata_oi->schan){
    return PKT_OBS_SCHAN;
  }
  if(pkt_schan + ata_oi->pkt_nchan > ata_oi->schan + ata_oi->nchan){
    return PKT_OBS_NCHAN;
  }
  return PKT_OBS_OK;
}

unsigned check_pkt_observability(
    const struct ata_snap_obs_info * ata_oi,
    const uint64_t pkt_idx,
    const uint64_t obs_start_pktidx,
    const uint16_t feng_id,
    const uint16_t pkt_schan
  );
unsigned check_pkt_observability_silent(
    const struct ata_snap_obs_info * ata_oi,
    const uint64_t pkt_idx,
    const uint64_t obs_start_pktidx,
    const uint16_t feng_id,
    const uint16_t pkt_schan
  );

uint32_t update_stt_status_keys( hashpipe_status_t *st,
                                    enum run_states state,
                                    uint64_t pktidx,
                                    struct mjd_t *mjd);

enum run_states state_from_start_stop(const uint64_t pktidx,
                                       const uint64_t obs_start_pktidx, const uint64_t obs_stop_pktidx);

enum run_states state_from_block_start_stop(const uint64_t obs_start_pktidx, const uint64_t obs_stop_pktidx,
                                       const uint64_t block_start_pktidx, const uint64_t block_stop_pktidx);

int ata_snap_obs_info_read(hashpipe_status_t *st, struct ata_snap_obs_info *obs_info);
int ata_snap_obs_info_write(hashpipe_status_t *st, struct ata_snap_obs_info *obs_info);
char ata_snap_obs_info_read_with_validity(hashpipe_status_t *st, struct ata_snap_obs_info *obs_info, enum obs_info_validity *validity);
void ata_snap_obs_info_write_with_validity(hashpipe_status_t *st, struct ata_snap_obs_info *obs_info, enum obs_info_validity obs_info_valid);

char align_blk0_with_obsstart(uint64_t * blk0_start_pktidx, uint32_t obsstart, uint32_t pktidx_per_block);

//
// The transposition is from a block of ATA SNAP packets,
// the headers of which specify:
//    PKTIDX[0 ... PIPERBLK]  (Packet-Time)
//    FENG  [0 ... NANT]      (AntennaEnum)
//
// where each SNAP packet has dimensions:
//    [Slowest ------> Fastest]
//    [PKTCHAN, PKTNTIME, NPOL]
// 
// to GUPPI RAW:
//    [Slowest ---> Fastest]
//    [Frequency, Time, Pol]
//
// The transposition takes each packet's [PKTNTIME, NPOL] slice
// and places it correctly in the overarching block's time
// and frequency dimension:
//
//    [FENG,    PKTNCHAN,   PKTIDX, PKTNTIME,   NPOl]
//    [Antenna, ^Frequency^, ^-----Time-----^, ^Pol^ ]
//
// This is costly because the largest slice that can be
// copied to its final position is [PKTNTIME, NPOL],
// which is [16, 2] = 32 bytes. This has to be repeated
// for each channel in each packet (PKTNCHAN * PIPERBLK)
//

// Packet content constants
#define ATASNAP_DEFAULT_SAMPLE_WIDTH_T uint16_t // this is the total width of the complex sample (8+8i = 16bit)
#define ATASNAP_DEFAULT_SAMPLE_BYTESIZE sizeof(ATASNAP_DEFAULT_SAMPLE_WIDTH_T)
#define ATASNAP_DEFAULT_PKT_SAMPLE_BYTE_STRIDE ATASNAP_DEFAULT_PKTNPOL*ATASNAP_DEFAULT_SAMPLE_BYTESIZE // this assumes that a packet's PKTIDX (ie timestamp) field increments in steps of NTIME
#define ATASNAP_DEFAULT_PKT_CHAN_BYTE_STRIDE ATASNAP_DEFAULT_PKTNTIME*ATASNAP_DEFAULT_PKT_SAMPLE_BYTE_STRIDE
typedef struct __attribute__ ((__packed__)) {ATASNAP_DEFAULT_SAMPLE_WIDTH_T num[ATASNAP_DEFAULT_PKTNPOL*ATASNAP_DEFAULT_PKTNTIME];} PKT_DCP_FTP_T; // sizeof(PKT_DCP_T) == ATASNAP_DEFAULT_PKT_CHAN_BYTE_STRIDE

#ifndef ATA_PACKET_PAYLOAD_DIRECT_COPY
#define COPY_PACKET_DATA_TO_FTP_DATABUF_FORLOOP(\
        /*const uint8_t**/  payload_dest,/*Indexed into [FENG, PKT_SCHAN, PKTIDX, 0, 0]*/\
        /*const uint8_t**/  pkt_payload,\
        /*const uint16_t*/  pkt_nchan,\
        /*const uint32_t*/  channel_stride, /*= PIPERBLK*ATASNAP_DEFAULT_PKTIDX_STRIDE/ATASNAP_DEFAULT_PKT_CHAN_BYTE_STRIDE */\
        /*const uint32_t*/  time_stride /* Unused as copy strides TIME*POLE */\
      )\
    for(int pkt_chan_idx = 0; pkt_chan_idx < pkt_nchan; pkt_chan_idx++){\
      memcpy(\
        payload_dest + channel_stride*pkt_chan_idx,\
        pkt_payload + pkt_chan_idx*ATASNAP_DEFAULT_PKT_CHAN_BYTE_STRIDE,\
        ATASNAP_DEFAULT_PKT_CHAN_BYTE_STRIDE\
      );\
    }
#else
#define COPY_PACKET_DATA_TO_FTP_DATABUF_FORLOOP(\
        /*PKT_DCP_T**/  payload_dest,/*Indexed into [FENG, PKT_SCHAN, PKTIDX, 0, 0]*/\
        /*PKT_DCP_T**/  pkt_payload,\
        /*const uint16_t*/  pkt_nchan,\
        /*const uint32_t*/  channel_stride, /*= PIPERBLK*ATASNAP_DEFAULT_PKTIDX_STRIDE/ATASNAP_DEFAULT_PKT_CHAN_BYTE_STRIDE */\
        /*const uint32_t*/  time_stride /* Unused as copy strides TIME*POLE */\
      )\
    for(int pkt_chan_idx = 0; pkt_chan_idx < pkt_nchan; pkt_chan_idx++){\
      *(payload_dest) = *pkt_payload++; \
      payload_dest += channel_stride;\
    }
#endif // define COPY_PACKET_DATA_TO_FTP_DATABUF_FORLOOP

typedef struct __attribute__ ((__packed__)) {ATASNAP_DEFAULT_SAMPLE_WIDTH_T num[ATASNAP_DEFAULT_PKTNPOL];} PKT_DCP_TFP_T;

// 
// to xGPU-Correlator input:
//    [Slowest ---> Fastest]
//    Time        [0 ... PIPERBLK*PKTNTIME]
//    Channel     [0 ... NCHAN]
//    FENG        [0 ... NANT]
//    POL         [0 ... NPOL]
//
// The transposition takes each NPOL pols together, i.e. 2x (8re+8im)
//
#ifndef ATA_PACKET_PAYLOAD_DIRECT_COPY
#define COPY_PACKET_DATA_TO_TFP_DATABUF_FORLOOP(\
        /*const uint8_t**/  payload_dest,/*Indexed into [PKTIDX, PKT_SCHAN, FENG, 0]*/\
        /*const uint8_t**/  pkt_payload,\
        /*const uint16_t*/  pkt_nchan,\
        /*const uint32_t*/  channel_stride,\
        /*const uint32_t*/  time_stride\
      )\
    for(int pkt_npol_sample_idx = 0; pkt_npol_sample_idx < pkt_nchan*ATASNAP_DEFAULT_PKTNTIME; pkt_npol_sample_idx++){ \
      memcpy(\
        payload_dest + \
            (pkt_npol_sample_idx/ATASNAP_DEFAULT_PKTNTIME)*channel_stride + (pkt_npol_sample_idx%ATASNAP_DEFAULT_PKTNTIME)*time_stride,\
        pkt_payload + pkt_npol_sample_idx*ATASNAP_DEFAULT_PKTNPOL*ATASNAP_DEFAULT_SAMPLE_BYTESIZE, \
        ATASNAP_DEFAULT_PKTNPOL*ATASNAP_DEFAULT_SAMPLE_BYTESIZE \
      );\
    }
#else
#define COPY_PACKET_DATA_TO_TFP_DATABUF_FORLOOP(\
        /*PKT_DCP_TFP_T**/  payload_dest,/*Indexed into [PKTIDX, PKT_SCHAN, FENG, 0]*/\
        /*PKT_DCP_TFP_T**/  pkt_payload,\
        /*const uint16_t*/  pkt_nchan,\
        /*const uint32_t*/  channel_stride,\
        /*const uint32_t*/  time_stride\
      )\
		for(int pkt_chan_idx = 0; pkt_chan_idx < pkt_nchan; pkt_chan_idx++){\
      for(int pkt_timeXnpol_idx = 0; pkt_timeXnpol_idx < ATASNAP_DEFAULT_PKTNTIME; pkt_timeXnpol_idx++){\
        *(payload_dest + pkt_chan_idx*channel_stride + pkt_timeXnpol_idx*time_stride) = *pkt_payload++;\
      }\
    }
#endif // define COPY_PACKET_DATA_TO_TFP_DATABUF_FORLOOP

// to xGPU(DP4A)-Correlator input:
//    [Slowest ---> Fastest]
//    Time        [0 ... PIPERBLK/4]
//    Channel     [0 ... NCHAN]
//    FENG        [0 ... NANT]
//    POL         [0 ... NPOL]
//    complexity  [real, imag]
//    time_minor  [0 ... 4]
//
// The transposition copies each byte, i.e. half of each sample (8re, 8im)

#ifdef __SSSE3__
#include <tmmintrin.h>
// https://stackoverflow.com/a/35268748
#define CHAR_AS_LONGLONG(a) (((long long)a) & 0xFF)
#define LL_SETR_EPI8(a, b, c, d, e, f, g, h) \
    CHAR_AS_LONGLONG(a) | (CHAR_AS_LONGLONG(b) << 8) | \
    (CHAR_AS_LONGLONG(c) << 16) | (CHAR_AS_LONGLONG(d) << 24) | \
    (CHAR_AS_LONGLONG(e) << 32) | (CHAR_AS_LONGLONG(f) << 40) | \
    (CHAR_AS_LONGLONG(g) << 48) | (CHAR_AS_LONGLONG(h) << 56)
#define _MM_SETR_EPI8(a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, aa, ab, ac, ad, ae, af) \
    {LL_SETR_EPI8(a0, a1, a2, a3, a4, a5, a6, a7), LL_SETR_EPI8(a8, a9, aa, ab, ac, ad, ae, af)}
// from
//    TIME        [0 ... 4]
//    POL         [0 ... NPOL]
//    complexity  [real, imag]
// C0T0P0R C0T0P0I C0T0P1R C0T0P1I			000:031
// C0T1P0R C0T1P0I C0T1P1R C0T1P1I			032:063
// C0T2P0R C0T2P0I C0T2P1R C0T2P1I			064:095
// C0T3P0R C0T3P0I C0T3P1R C0T3P1I			096:127
// 	00 01 02 03
// 	04 05 06 07
// 	08 09 10 11
// 	12 13 14 15
//
// to DP4A
//    POL         [0 ... NPOL]
//    complexity  [real, imag]
//    time_minor  [0 ... 4]
// C0T0P0R C0T1P0R C0T2P0R C0T3P0R			000:031
// C0T0P0I C0T1P0I C0T2P0I C0T3P0I			032:063
// C0T0P1R C0T1P1R C0T2P1R C0T3P1R			064:095
// C0T0P1I C0T1P1I C0T2P1I C0T3P1I			096:127
// 	00 04 08 12
// 	01 05 09 13
// 	02 06 10 14
// 	03 07 11 15
static const __m128i CORNER_TURN_SHUFFLE_MASK = _MM_SETR_EPI8(
	0, 4,  8, 12,
	1, 5,  9, 13,
	2, 6, 10, 14,
	3, 7, 11, 15
);
typedef __m128i PKT_DCP_TFP_DP4A_T;

#define COPY_PACKET_DATA_TO_TFP_DP4A_DATABUF_FORLOOP(\
        /*PKT_DCP_TFP_DP4A_T**/  payload_dest,/*Indexed into [PKTIDX, PKT_SCHAN, FENG, 0]*/\
        /*PKT_DCP_TFP_DP4A_T**/  pkt_payload,\
        /*const uint16_t*/  pkt_nchan,\
        /*const uint32_t*/  channel_stride,\
        /*const uint32_t*/  time_stride\
      )\
    for(int pkt_chan_idx = 0; pkt_chan_idx < pkt_nchan; pkt_chan_idx++){\
      for(int pkt_time_major_idx = 0; pkt_time_major_idx < ATASNAP_DEFAULT_PKTNTIME/4; pkt_time_major_idx++){\
        payload_dest[\
          (pkt_time_major_idx * time_stride*4) +\
          pkt_chan_idx * channel_stride\
        ] = _mm_shuffle_epi8(*pkt_payload++, CORNER_TURN_SHUFFLE_MASK);\
      }\
    }
#else

#if ATA_PAYLOAD_TRANSPOSE == ATA_PAYLOAD_TRANSPOSE_TFP_DP4A
#pragma message("Enable SSSE3 for far better performance of TFP_DP4A unpack!") 
#endif

typedef uint8_t PKT_DCP_TFP_DP4A_T; // this is the width of the one component of the complex sample (8+8i = 16bit)/2 = 8bit

#define COPY_PACKET_DATA_TO_TFP_DP4A_DATABUF_FORLOOP(\
        /*PKT_DCP_TFP_DP4A_T**/  payload_dest,/*Indexed into [PKTIDX, PKT_SCHAN, FENG, 0]*/\
        /*PKT_DCP_TFP_DP4A_T**/  pkt_payload,\
        /*const uint16_t*/  pkt_nchan,\
        /*const uint32_t*/  channel_stride,\
        /*const uint32_t*/  time_stride\
      )\
    for(int pkt_chan_idx = 0; pkt_chan_idx < pkt_nchan; pkt_chan_idx++){\
      for(int pkt_time_major_idx = 0; pkt_time_major_idx < ATASNAP_DEFAULT_PKTNTIME/4; pkt_time_major_idx++){\
        for(int pkt_time_minor_idx = 0; pkt_time_minor_idx < 4; pkt_time_minor_idx++){\
          for(int pkt_pol_idx = 0; pkt_pol_idx < ATASNAP_DEFAULT_PKTNPOL; pkt_pol_idx++){\
            for(int c = 0; c < 2; c++){\
              payload_dest[\
                (pkt_time_major_idx * time_stride*4) +\
                pkt_chan_idx * channel_stride +\
                pkt_pol_idx*4*ATASNAP_DEFAULT_SAMPLE_BYTESIZE +\
                c*4*ATASNAP_DEFAULT_SAMPLE_BYTESIZE/2 +\
                pkt_time_minor_idx*ATASNAP_DEFAULT_SAMPLE_BYTESIZE/2\
              ] = *pkt_payload++;\
            }\
          }\
        }\
      }\
    }
#endif // define COPY_PACKET_DATA_TO_TFP_DATABUF_FORLOOP

#endif // _HPGUPPI_ATASNAP_H_
