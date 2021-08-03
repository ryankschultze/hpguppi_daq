#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#include <omp.h>

#include "hashpipe.h"
#include "hpguppi_databuf.h"

// Packet content constants
#define NPOL  2
#define NTIME 16
#define NBITS 4
#define NT_MINOR 128/NBITS
#define prod_NPOL_NT_MINOR NPOL * NT_MINOR 

#define TIME_WIDTH_T uint8_t // this is the total width of the complex sample (4+4i = 8bit)
typedef struct {TIME_WIDTH_T num[NPOL*NTIME];} CP_DTYPE;

typedef uint16_t CP_DTYPE_XGPU; //This is for xGPU

#define USE_MULTI_THREAD
#define MULTI_THREAD_COUNT 8

typedef struct
{
  // int ntime; // Static per hashpipe for memcpy performance gains
  int obsnchan;
  int nbits;
  int ndim;
  // int npol; // Static per hashpipe for memcpy performance gains
  int piperblk;
  int pktnchan; //XXX
  int nants; //XXX
  int nstrm; //XXX
} db_transpose_t;

//
// The transposition is from a block of contiguous ATA SNAP packets,
// laid out with the order:
//    PKTIDX[0 ... PIPERBLK]  (Packet-Time)  Slowest
//    FENG  [0 ... NANT]      (AntennaEnum)  
//    STREAM[0 ... NSTRM]     (Packet-Freq)  Fastest
//
// where each SNAP packet has dimensions:
//    [Slowest ------> Fastest]
//    [PKTCHAN, PKTNTIME, NPOL]
// 
// to xGPU-Correlator input:
//    [Slowest ---> Fastest]
//    Time        [0 ... PIPERBLK*PKTNTIME]
//    Channel     [0 ... NSTRM*PKTNCHAN]
//    FENG        [0 ... NANT]
//    POL         [0 ... NPOL]
//
// The transposition takes each NPOL pols together, i.e. 2x (4re+4im)
// This is treated as a unint16_t type with CP_DTYPE_XGPU
//
int transpose_to_xgpu(db_transpose_t * ctx, const void* in, void* out_void)
{
  // number of packets that span the entire data block, in time
  size_t itime_packets = ctx->piperblk / NTIME;
  size_t pktchs = ctx->pktnchan;
  // size_t pktchs = ctx->obsnchan / ctx->nants;
  //size_t nchan = ctx->obsnchan;
  //
  size_t pktidx, feng, strm, pkt_ch, pkt_tm;
  size_t t,ch;
  size_t ostride=0;

  CP_DTYPE_XGPU *inp = (CP_DTYPE_XGPU*) in;
  CP_DTYPE_XGPU *out = (CP_DTYPE_XGPU*) out_void;
  CP_DTYPE_XGPU *inp_p;

  inp_p = inp;


#ifdef USE_MULTI_THREAD
#pragma omp parallel for private (inp_p)
#endif
  for (pktidx=0; pktidx<itime_packets; pktidx++)
  {
    inp_p = inp + pktidx*ctx->nants*ctx->nstrm*pktchs*NTIME;
    for (feng=0; feng<ctx->nants; feng++)
    {
      for (strm=0; strm<ctx->nstrm; strm++)
      {
        for (pkt_ch=0; pkt_ch<pktchs; pkt_ch++)
        {
          ch = strm*pktchs + pkt_ch;
          for (pkt_tm=0; pkt_tm<NTIME; pkt_tm++)
          {
            t = pktidx * NTIME + pkt_tm;
            ostride = t*(ctx->nstrm*pktchs)*ctx->nants
              + ch*ctx->nants + feng;

            out[ostride] = *inp_p++;
            // The following is to check whether input is properly
            // mapped to output
            //*inp_p = out[ostride];
            //inp_p++;
          }
        }
      }
    }
  }

  return 0;
}

static void *run(hashpipe_thread_args_t *args)
{
  hpguppi_input_databuf_t *indb  = (hpguppi_input_databuf_t *)args->ibuf;
  hpguppi_input_databuf_t *outdb = (hpguppi_input_databuf_t *)args->obuf;

  hashpipe_status_t st = args->st;
  const char* status_key = args->thread_desc->skey;
  const char* thread_name = args->thread_desc->name;
  char* databuf_header;

  int rv;

  // int correlate_not_transpose_switch=0;
  int curblock_in=0;
  int curblock_out=0;

  db_transpose_t ctx;

  
#ifdef USE_MULTI_THREAD
    // Each thread processes about 11.62 Gbits/s
    const int nthreads = MULTI_THREAD_COUNT;// Reach for 32  Gbits/s

    hashpipe_status_lock_safe(&st);
      hputi4(st.buf, "TRNTHRDS", nthreads);
    hashpipe_status_unlock_safe(&st);
    omp_set_num_threads(nthreads);
#endif

  while (run_threads())
  {
    hashpipe_status_lock_safe(&st);
      hputi4(st.buf, "TRBLKIN", curblock_in);
      hputs(st.buf, status_key, "waiting");
      hputi4(st.buf, "TRBLKOUT", curblock_out);
    hashpipe_status_unlock_safe(&st);

    // Waiting for input
    while ((rv=hpguppi_input_databuf_wait_filled(indb, curblock_in)) != 
		    HASHPIPE_OK)
    {
      if (rv == HASHPIPE_TIMEOUT)
      {
        hashpipe_status_lock_safe(&st);
	      hputs(st.buf, status_key, "inblocked");
    	  hashpipe_status_unlock_safe(&st);
	      continue;
      }
      else
      {
        hashpipe_error(thread_name, "error waiting for input buffer, rv: %i", rv);
	      //pthread_exit(NULL);
	      //break;
	continue;
      }
    }

    // Waiting for output
    while ((rv=hpguppi_input_databuf_wait_free(outdb, curblock_out)) !=
		    HASHPIPE_OK)
    {
      if (rv == HASHPIPE_TIMEOUT)
      {
        hashpipe_status_lock_safe(&st);
	      hputs(st.buf, status_key, "outblocked");
       	hashpipe_status_unlock_safe(&st);
	      continue;
      }
      else
      {
        hashpipe_error(thread_name, "error waiting for output buffer, rv: %i", rv);
        //pthread_exit(NULL);
	      //break;
	      continue;
      }
    }

    hashpipe_status_lock_safe(&st);
    hputs(st.buf, status_key, "transposing");
    hashpipe_status_unlock_safe(&st);

    // create context
    databuf_header = hpguppi_databuf_header(indb, curblock_in);
    hgeti4(databuf_header, "OBSNCHAN", &ctx.obsnchan);
    hgeti4(databuf_header, "PIPERBLK", &ctx.piperblk);
    hgeti4(databuf_header, "PKTNCHAN", &ctx.pktnchan);
    hgeti4(databuf_header, "NANTS", &ctx.nants);
    hgeti4(databuf_header, "NSTRM", &ctx.nstrm);
  
    // copy across the header
    memcpy(hpguppi_databuf_header(outdb, curblock_out), 
           databuf_header, 
           BLOCK_HDR_SIZE);//HASHPIPE_STATUS_TOTAL_SIZE);

    // hgeti4(databuf_header, "XSWITCH", &correlate_not_transpose_switch);
    // switch(correlate_not_transpose_switch){
    //   case 1:
        transpose_to_xgpu(&ctx, hpguppi_databuf_data(indb, curblock_in), 
            hpguppi_databuf_data(outdb, curblock_out));
    //     break;
    //   default:
    //     transpose(&ctx, hpguppi_databuf_data(indb, curblock_in), 
    //         hpguppi_databuf_data(outdb, curblock_out));
    //     break;
    // }

    hpguppi_input_databuf_set_free(indb, curblock_in);
    curblock_in  = (curblock_in + 1) % indb->header.n_block;

    hpguppi_input_databuf_set_filled(outdb, curblock_out);
    curblock_out = (curblock_out + 1) % outdb->header.n_block;
  }

  hashpipe_info(thread_name, "returning");

  return NULL;
}

static hashpipe_thread_desc_t pkt_xgpu_transpose_thread = {
  name: "hpguppi_atasnap_pkt_to_xgpu_transpose",
  skey: "TRANSTAT",
  init: NULL,
  run: run,
  ibuf_desc: {hpguppi_input_databuf_create},
  obuf_desc: {hpguppi_input_databuf_create}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(& pkt_xgpu_transpose_thread);
}
