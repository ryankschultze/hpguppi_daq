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
#define TIME_WIDTH_T uint8_t // this is the total width of the complex sample (4+4i = 8bit)
typedef struct {TIME_WIDTH_T num[NPOL*NTIME];} CP_DTYPE;

#define USE_MULTI_THREAD

typedef struct
{
  int ntime;
  int obsnchan;
  int nbits;
  int ndim;
  int npol;
  int piperblk;
} db_transpose_t;

int transpose(db_transpose_t * ctx, const void* in, void* out)
{
  // To be used in pointer arithmetic
  const CP_DTYPE* inbuf;
  CP_DTYPE* outbuf;

  const CP_DTYPE* baseinbuf = in;
  CP_DTYPE* baseoutbuf      = out;

  // number of packets that span the entire data block, in time
  size_t itime_packets = ctx->piperblk / ctx->ntime;
  size_t nchan = ctx->obsnchan;

  // bytes to stride within a packet in input buffer
  // and also amount to copy at a time
  //size_t istride = (ctx->npol * ctx->ndim * ctx->nbits * ctx->ntime)/8/sizeof *inbuf;
  size_t istride = NPOL*NTIME/sizeof *inbuf;//istride should work out to 1 with a correct CP_DTYPE
  //size_t istride_m = istride*sizeof *inbuf;

#ifdef USE_MULTI_THREAD
  size_t tstride = ctx->obsnchan * istride;
#endif
  size_t ostride = itime_packets * istride;
  size_t nstride = ostride - istride;

  // Loop over entire spectrum-packets over all the chans
#ifdef USE_MULTI_THREAD
#pragma omp parallel for private (inbuf, outbuf)
#else
  inbuf  = baseinbuf;// + iptime*tstride;
#endif
  for (size_t iptime=0; iptime < itime_packets; iptime++)
  {
#ifdef USE_MULTI_THREAD
    inbuf  = baseinbuf + iptime*tstride;
#endif
    outbuf = baseoutbuf + iptime*istride;
    for (size_t ichan=0; ichan < nchan; ichan+=1)
    {
      //memcpy(outbuf, inbuf, istride_m);
      //inbuf += istride;
      //outbuf += ostride;
      for (size_t t=0; t<istride; t++)
        *outbuf++ = *inbuf++;// benchmarks show this to be faster than memcpy
      outbuf += nstride;
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

  int rv;

  int curblock_in=0;
  int curblock_out=0;

  db_transpose_t ctx;

  
#ifdef USE_MULTI_THREAD
    // Each thread processes about 11.62 Gbits/s
    const int nthreads = 3;// Reach for 32  Gbits/s

    hputi4(st.buf, "TRNTHRDS", nthreads);
    omp_set_num_threads(nthreads);
#endif

  while (run_threads())
  {
    hashpipe_status_lock_safe(&st);
    hputi4(st.buf, "TRBLKIN", curblock_in);
    hputs(st.buf, status_key, "waiting");
    hputi4(st.buf, "TRBLKOUT", curblock_out);

    hgeti4(st.buf, "OBSNCHAN", &ctx.obsnchan);
    hgeti4(st.buf, "PIPERBLK", &ctx.piperblk);
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
        hashpipe_error(thread_name, "error waiting for input buffer");
	      pthread_exit(NULL);
	      break;
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
        hashpipe_error(thread_name, "error waiting for output buffer");
        pthread_exit(NULL);
	      break;
      }
        
    }

    hashpipe_status_lock_safe(&st);
    hputs(st.buf, status_key, "transposing");

    hashpipe_status_unlock_safe(&st);
    // create context
  
    // copy across the header
    memcpy(hpguppi_databuf_header(outdb, curblock_out), 
           hpguppi_databuf_header(indb, curblock_in), 
           HASHPIPE_STATUS_TOTAL_SIZE);

    transpose(&ctx, hpguppi_databuf_data(indb, curblock_in), 
		    hpguppi_databuf_data(outdb, curblock_out));
    
    hpguppi_input_databuf_set_free(indb, curblock_in);
    curblock_in  = (curblock_in + 1) % indb->header.n_block;

    hpguppi_input_databuf_set_filled(outdb, curblock_out);
    curblock_out = (curblock_out + 1) % outdb->header.n_block;
  }

  hashpipe_info(thread_name, "returning");

  return NULL;
}





static hashpipe_thread_desc_t pkt_FTP_transpose_thread = {
  name: "hpguppi_atasnap_pkt_to_FTP_transpose",
  skey: "TRANSTAT",
  init: NULL,
  run: run,
  ibuf_desc: {hpguppi_input_databuf_create},
  obuf_desc: {hpguppi_input_databuf_create}
};

static __attribute__((constructor)) void ctor()
{
  register_hashpipe_thread(& pkt_FTP_transpose_thread);
}
