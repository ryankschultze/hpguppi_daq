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

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#include "hashpipe.h"
#include "hpguppi_databuf.h"
#include "hpguppi_time.h"
#include "hpguppi_util.h"
#include "hashpipe_packet.h"

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
                    PKT_OBS_STREAM
                  };
enum obs_info_validity { 
                OBS_UNKNOWN=-5,       //=-5
                OBS_INVALID_FENG=-4,  //=-4
                OBS_INVALID_SCHAN=-3, //=-3
                OBS_INVALID_STREAM=-2,//=-2
                OBS_INVALID=-1,       //=-1
                OBS_SEEMS_VALID=0,    //=0   sign bit for validity
                OBS_VALID=1           //=1
              };

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

// ATA SNAP packet with link layer header and internal padding to optimize
// alignment.  The alignment is acheived through judicious use of IB Verbs
// scatter/gather capabilities (specifically the scatter part).
struct __attribute__ ((__packed__)) ata_snap_ibv_pkt {
  struct ethhdr ethhdr; //0
  struct iphdr iphdr;   //14
  struct udphdr udphdr; //34
  uint8_t pad0[22];     //42
  uint8_t version;      //64
  uint8_t type;         //65
  uint16_t n_chans;     //66
  uint16_t chan;        //68
  uint16_t feng_id;     //70
  uint64_t timestamp;   //72
  uint8_t pad1[48];     //80
	uint8_t payload[];    //128
};

struct __attribute__ ((__packed__)) ata_snap_pkt {
  struct ethhdr ethhdr; // headers aren't functional  // 14 bytes
  struct iphdr iphdr; // headers aren't functional    // 20 bytes
  struct udphdr udphdr; // headers aren't functional  // 8 bytes
  uint8_t pad0[66]; // Dont know the reason, but actual network payload is offset
  uint8_t version;
  uint8_t type;
  uint16_t n_chans;
  uint16_t chan;
  uint16_t feng_id;
  uint64_t timestamp;
  uint8_t payload[];//complex4 data[n_chans, 16, 2] // 4-bit real + 4-bit imaginary
};

#define ATA_SNAP_PKT_NUMBER(ata_snap_pkt)   (uint64_t)__bswap_64(ata_snap_pkt->timestamp)
#define ATA_SNAP_PKT_CHAN(ata_snap_pkt)     __bswap_16(ata_snap_pkt->chan)
#define ATA_SNAP_PKT_FENG_ID(ata_snap_pkt)  __bswap_16(ata_snap_pkt->feng_id)

// ATA SNAP header byte offset within (unpadded) packet
#define ATA_SNAP_PKT_OFFSET_HEADER \
  (sizeof(struct ethhdr) + \
   sizeof(struct iphdr ) + \
   sizeof(struct udphdr))

#define ATA_SNAP_PKT_SIZE_HEADER (2*sizeof(uint64_t))

// ATA SNAP payload byte offset within (unpadded) packet
#define ATA_SNAP_PKT_OFFSET_PAYLOAD \
  (ATA_SNAP_PKT_OFFSET_HEADER + ATA_SNAP_PKT_SIZE_HEADER)

#define ATA_SNAP_PKT_SIZE_PAYLOAD (8192)

// The ATA SNAP design send out 8192 bytes of payload per packet.  MAX_PKT_SIZE
// is the next power of two that fits that plus non-payload bytes.
// TODO Why next power of two???
#define MAX_PKT_SIZE (16384)

struct ata_snap_feng_info {
  uint64_t pktidx;
  uint64_t feng_id;
  uint64_t feng_chan;
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
  // Total number of streams being processed
  uint32_t nstrm;
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
#define ATASNAP_DEFAULT_NSTRM           (     1)
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
  poi->nstrm = ATASNAP_DEFAULT_NSTRM;
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
    (oi.nstrm      != 0) &&
    (oi.pkt_npol   != 0) &&
    (oi.time_nbits != 0) &&
    (oi.pkt_ntime  != 0) &&
    (oi.pkt_nchan  != 0) &&
    (oi.schan      != OBS_INFO_INVALID_SCHAN);
}

// For ATA SNAP, the OBSNCHAN parameter represents the total number of
// frequency channels processed.  It is the number of antennas times the number
// of streams per antenna times the number of channels per packet/stream.
static inline
uint32_t
calc_ata_snap_obsnchan(uint32_t nants, uint32_t nstrm, uint32_t pkt_nchan)
{
  return nants * nstrm * pkt_nchan;
}

static inline
uint32_t
ata_snap_obsnchan(const struct ata_snap_obs_info oi)
{
  return calc_ata_snap_obsnchan(oi.nants, oi.nstrm, oi.pkt_nchan);
}

// This function doesn't assume each sample is 2 bytes ([4 bits real + 4 bits imag] *
// 2 pols).
static inline
uint32_t
calc_ata_snap_pkt_payload_bytes(uint32_t nchan,
                               uint32_t pkt_ntime, uint32_t pkt_npol,
                               uint32_t time_nbits)
{ // each time sample is real
  return (nchan * pkt_ntime * pkt_npol * 2 * time_nbits) / 8;
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
//
// Furthermore, the number of time samples per block is desired to be a power
// of two, so that subsequent FFTs can operate on complete blocks with maximum
// efficiency.  For this to happen, both the number of time samples per packet
// (PKT_NTIME) and the number of PKTIDX values per block (PIPERBLK) must be
// powers of two.  PKT_NTIME is set by the upstream F engines, so we have no
// control of that, but we can and do ensure that PIPERBLK is a power of 2.
static inline
uint32_t
calc_ata_snap_pkt_per_block(size_t block_size, uint32_t nchan,
                               uint32_t pkt_ntime, uint32_t pkt_npol,
                               uint32_t time_nbits)
{
  uint32_t pkt_idxpblk = (uint32_t)(block_size /
      calc_ata_snap_pkt_payload_bytes(nchan, pkt_ntime, pkt_npol, time_nbits));
  return prevpow2(pkt_idxpblk);
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
calc_ata_snap_eff_pkt_per_block(size_t block_size, uint32_t nchan,
                               uint32_t pkt_ntime, uint32_t pkt_npol,
                               uint32_t time_nbits, uint32_t nants,
                               uint32_t nstrm)
{
  // use integer division to round down to nearest multiple of (nants*nstrm)
  uint32_t total_strms = nants * nstrm;
  return total_strms *
        (calc_ata_snap_pkt_per_block(block_size, nchan, pkt_ntime, pkt_npol, time_nbits)
        / total_strms);
}

static inline
uint32_t
ata_snap_eff_pkt_per_block(size_t block_size, const struct ata_snap_obs_info oi)
{
  return calc_ata_snap_eff_pkt_per_block(
                  block_size, oi.pkt_nchan, oi.pkt_ntime,
                  oi.pkt_npol, oi.time_nbits,
                  oi.nants, oi.nstrm);
}

static inline
uint32_t
calc_ata_snap_pkt_bytes(uint32_t nchan,
                        uint32_t pkt_ntime, uint32_t pkt_npol,
                        uint32_t time_nbits)
{
  return calc_ata_snap_pkt_payload_bytes(nchan, pkt_ntime, pkt_npol, time_nbits) + 16;
}

static inline
uint32_t
ata_snap_pkt_bytes(const struct ata_snap_obs_info oi)
{ 
  return calc_ata_snap_pkt_bytes(oi.pkt_nchan, oi.pkt_ntime, oi.pkt_npol, oi.time_nbits);
}

// For ATA SNAP, the NTIME parameter (not stored in the status buffer or GUPPI
// RAW files), represents the total number of time samples in a block.  It
// depends on the block size and NCHAN (and sample size, but that is assumed to
// be 2 bytes).  This calculation is a bit tricky because the effective block
// size can be less than the max block size when NCHAN and PKT_NTIME do not
// evenly divide the max block size.
static inline
uint32_t
calc_ata_snap_pktidx_per_block(size_t block_size, uint32_t nchan, uint32_t pkt_ntime,
                            uint32_t pkt_npol, uint32_t time_nbits,
                            uint32_t nants, uint32_t nstrm)
{
  uint32_t pkidx_per_block = (calc_ata_snap_eff_pkt_per_block(block_size, 
                                            nchan, pkt_ntime, pkt_npol, time_nbits,
                                            nants, nstrm)
                              /(nants*nstrm)) * pkt_ntime;
  return pkidx_per_block;
}

static inline
uint32_t
ata_snap_pktidx_per_block(size_t block_size, const struct ata_snap_obs_info oi)
{
  return calc_ata_snap_pktidx_per_block(block_size, oi.pkt_nchan,
                                     oi.pkt_ntime, oi.pkt_npol, oi.time_nbits,
                                     oi.nants, oi.nstrm);
}

// Calculate the effective block size for the given max block size, nchan, and
// pkt_ntime values.  The effective block size can be less than the max block size
// if nchan and/or pkt_ntime do not evenly divide the max block size.
static inline
uint32_t
calc_ata_snap_block_size(size_t block_size, uint32_t pkt_nchan, uint32_t pkt_ntime,
                         uint32_t pkt_npol, uint32_t time_nbits, uint32_t nants, uint32_t nstrm)
{
  return    calc_ata_snap_pkt_payload_bytes(pkt_nchan, pkt_ntime,
                                            pkt_npol, time_nbits)
            * calc_ata_snap_eff_pkt_per_block(
                  block_size, pkt_nchan, pkt_ntime,
                  pkt_npol, time_nbits,
                  nants, nstrm);
}

static inline
uint32_t
ata_snap_block_size(size_t block_size, const struct ata_snap_obs_info oi)
{
  return calc_ata_snap_block_size(block_size, oi.pkt_nchan,
                                     oi.pkt_ntime, oi.pkt_npol, oi.time_nbits,
                                     oi.nants, oi.nstrm);
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

static inline
void
ata_snap_populate_block_related_fields(size_t block_size, struct ata_snap_obs_info *oi)
{
  oi->pkt_data_size = ata_snap_pkt_bytes(*oi);
  oi->pkt_per_block = ata_snap_eff_pkt_per_block(block_size, *oi);
  oi->pktidx_per_block = ata_snap_pktidx_per_block(block_size, *oi);//inherently effective 
  oi->eff_block_size = oi->pkt_per_block*(oi->pkt_data_size-16);
}

// Parses a ATA SNAP F Engine packet that is in the format of a "struct
// ata_snap_ibv_pkt", stores metadata in fei, returns pointer to packet's
// payload.
static inline
const uint8_t *
ata_snap_parse_ibv_packet(const struct ata_snap_ibv_pkt *p,
    struct ata_snap_feng_info * fei)
{
  fei->pktidx = ATA_SNAP_PKT_NUMBER(p);
  fei->feng_id = ATA_SNAP_PKT_FENG_ID(p);
  fei->feng_chan = ATA_SNAP_PKT_CHAN(p);

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


char * datablock_stats_data(const struct datablock_stats *d);
char * datablock_stats_header(const struct datablock_stats *d);
void reset_datablock_stats(struct datablock_stats *d);
void init_datablock_stats(struct datablock_stats *d,
    struct hpguppi_input_databuf *dbout, int block_idx, int64_t block_num,
    uint64_t pkts_per_block);
void block_stack_push(struct datablock_stats *d, int nblock);
void finalize_block(struct datablock_stats *d);
void increment_block(struct datablock_stats *d, int64_t block_num);
void wait_for_block_free(const struct datablock_stats * d,
  hashpipe_status_t * st, const char * status_key);

unsigned check_pkt_observability_sans_idx(
    const struct ata_snap_obs_info * ata_oi,
    const uint16_t feng_id,
    const int32_t stream,
    const uint16_t pkt_schan
  );
unsigned check_pkt_observability(
    const struct ata_snap_obs_info * ata_oi,
    const uint64_t pkt_idx,
    const uint64_t obs_start_pktidx,
    const uint16_t feng_id,
    const int32_t stream,
    const uint16_t pkt_schan
  );
unsigned check_pkt_observability_silent(
    const struct ata_snap_obs_info * ata_oi,
    const uint64_t pkt_idx,
    const uint64_t obs_start_pktidx,
    const uint16_t feng_id,
    const int32_t stream,
    const uint16_t pkt_schan
  );

void update_stt_status_keys( hashpipe_status_t *st,
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
#if 0
#define COPY_PACKET_DATA_TO_DATABUF(\
      /*const struct datablock_stats*/    datablock_stats_pointer,\
      /*const uint8_t*/   pkt_payload,\
      /*const uint64_t*/  pkt_obs_relative_idx,\
      /*const uint16_t*/  feng_id,\
      /*const int32_t*/   stream,\
      /*const uint16_t*/  pkt_schan,\
      /*const uint32_t*/  fid_stride,\
      /*const uint32_t*/  time_stride,\
      /*const uint64_t*/  pkt_payload_size,\
      /*const uint32_t*/  pkt_ntime)\
  memcpy(datablock_stats_data(datablock_stats_pointer)+(\
        (pkt_obs_relative_idx/pkt_ntime) * time_stride\
            +  feng_id * fid_stride\
            +  stream * pkt_payload_size\
      ),\
      pkt_payload, 4096)
#else
  #define COPY_PACKET_DATA_TO_DATABUF_1024(\
        /*const struct datablock_stats*/    datablock_stats_pointer,\
        /*const uint8_t*/   pkt_payload,\
        /*const uint64_t*/  pkt_obs_relative_idx,\
        /*const uint16_t*/  feng_id,\
        /*const int32_t*/   stream,\
        /*const uint16_t*/  pkt_schan,\
        /*const uint32_t*/  fid_stride,\
        /*const uint32_t*/  time_stride,\
        /*const uint32_t*/  pkt_ntime)\
    memcpy(datablock_stats_data(datablock_stats_pointer)+(\
          (pkt_obs_relative_idx/pkt_ntime) * time_stride\
              +  feng_id * fid_stride\
              +  stream * 1024\
        ),\
        pkt_payload, 1024)
  #define COPY_PACKET_DATA_TO_DATABUF_2048(\
        /*const struct datablock_stats*/    datablock_stats_pointer,\
        /*const uint8_t*/   pkt_payload,\
        /*const uint64_t*/  pkt_obs_relative_idx,\
        /*const uint16_t*/  feng_id,\
        /*const int32_t*/   stream,\
        /*const uint16_t*/  pkt_schan,\
        /*const uint32_t*/  fid_stride,\
        /*const uint32_t*/  time_stride,\
        /*const uint32_t*/  pkt_ntime)\
    memcpy(datablock_stats_data(datablock_stats_pointer)+(\
          (pkt_obs_relative_idx/pkt_ntime) * time_stride\
              +  feng_id * fid_stride\
              +  stream * 2048\
        ),\
        pkt_payload, 2048)
  #define COPY_PACKET_DATA_TO_DATABUF_4096(\
        /*const struct datablock_stats*/    datablock_stats_pointer,\
        /*const uint8_t*/   pkt_payload,\
        /*const uint64_t*/  pkt_obs_relative_idx,\
        /*const uint16_t*/  feng_id,\
        /*const int32_t*/   stream,\
        /*const uint16_t*/  pkt_schan,\
        /*const uint32_t*/  fid_stride,\
        /*const uint32_t*/  time_stride,\
        /*const uint32_t*/  pkt_ntime)\
    memcpy(datablock_stats_data(datablock_stats_pointer)+(\
          (pkt_obs_relative_idx/pkt_ntime) * time_stride\
              +  feng_id * fid_stride\
              +  stream * 4096\
        ),\
        pkt_payload, 4096)
  #define COPY_PACKET_DATA_TO_DATABUF_8192(\
        /*const struct datablock_stats*/    datablock_stats_pointer,\
        /*const uint8_t*/   pkt_payload,\
        /*const uint64_t*/  pkt_obs_relative_idx,\
        /*const uint16_t*/  feng_id,\
        /*const int32_t*/   stream,\
        /*const uint16_t*/  pkt_schan,\
        /*const uint32_t*/  fid_stride,\
        /*const uint32_t*/  time_stride,\
        /*const uint32_t*/  pkt_ntime)\
    memcpy(datablock_stats_data(datablock_stats_pointer)+(\
          (pkt_obs_relative_idx/pkt_ntime) * time_stride\
              +  feng_id * fid_stride\
              +  stream * 8192\
        ),\
        pkt_payload, 8192)

  static inline int COPY_PACKET_DATA_TO_DATABUF(
        const struct datablock_stats*   datablock_stats_pointer,
        const uint8_t*  pkt_payload,
        const uint64_t  pkt_obs_relative_idx,
        const uint16_t  feng_id,
        const int32_t   stream,
        const uint16_t  pkt_schan,
        const uint32_t  fid_stride,
        const uint32_t  time_stride,
        const uint64_t  pkt_payload_size,
        const uint32_t  pkt_ntime)
  {
    switch(pkt_payload_size){
      case 1024:
        COPY_PACKET_DATA_TO_DATABUF_1024(datablock_stats_pointer,
          pkt_payload, pkt_obs_relative_idx, feng_id, stream,
          pkt_schan, fid_stride, time_stride, pkt_ntime
        );
        break;
      case 2048:
        COPY_PACKET_DATA_TO_DATABUF_2048(datablock_stats_pointer,
          pkt_payload, pkt_obs_relative_idx, feng_id, stream,
          pkt_schan, fid_stride, time_stride, pkt_ntime
        );
        break;
      case 4096:
        COPY_PACKET_DATA_TO_DATABUF_4096(datablock_stats_pointer,
          pkt_payload, pkt_obs_relative_idx, feng_id, stream,
          pkt_schan, fid_stride, time_stride, pkt_ntime
        );
        break;
      case 8192:
        COPY_PACKET_DATA_TO_DATABUF_8192(datablock_stats_pointer,
          pkt_payload, pkt_obs_relative_idx, feng_id, stream,
          pkt_schan, fid_stride, time_stride, pkt_ntime
        );
        break;
      default:
        memcpy(datablock_stats_data(datablock_stats_pointer)+(
          (pkt_obs_relative_idx/pkt_ntime) * time_stride
              +  feng_id * fid_stride
              +  stream * pkt_payload_size
        ),  pkt_payload, pkt_payload_size);
        hashpipe_warn(__FUNCTION__, "Packet_payload_size of %d is not optimally handled.", pkt_payload_size);
        return -1;
        break;
    }
    return 0;
  }
#endif // COPY_PACKET_DATA_TO_DATABUF
#endif // _HPGUPPI_ATASNAP_H_
