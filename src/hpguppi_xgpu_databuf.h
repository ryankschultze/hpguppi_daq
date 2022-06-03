/* hpguppi_databuf.h
 *
 * Defines shared mem structure for data passing.
 */
#ifndef _HPGUPPI_XGPU_DATABUF_H
#define _HPGUPPI_XGPU_DATABUF_H
#include "xgpu.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/sem.h>
#include <errno.h>
#include <time.h>

#include <hashpipe.h>

#include "hpguppi_databuf.h"

#define XGPU_BLOCK_DATA_SIZE 2*BLOCK_DATA_SIZE // in bytes, from guppi_daq_server

// #define XGPU_INTEGRATE_AS_CF64_ON_CPU

#ifdef XGPU_INTEGRATE_AS_CF64_ON_CPU
#define XGPU_OUTPUT_BLOCK_ELEMENT_T UVH5_CF64_t
#else
#define XGPU_OUTPUT_BLOCK_ELEMENT_T Complex
#endif

typedef struct hpguppi_input_xgpu_block {
  char hdr[BLOCK_HDR_SIZE];
  char data[XGPU_BLOCK_DATA_SIZE];
} hpguppi_input_xgpu_block_t;

typedef struct hpguppi_input_xgpu_databuf {
  hashpipe_databuf_t header;
  hashpipe_databuf_alignment padding; // Maintain data alignment
  hpguppi_input_xgpu_block_t block[N_INPUT_BLOCKS];
} hpguppi_input_xgpu_databuf_t;

hashpipe_databuf_t *hpguppi_input_xgpu_databuf_create(int instance_id, int databuf_id);

#define N_XGPU_OUTPUT_BLOCKS 16

//! this is dynamically sized
//! - don't access directly, use `hpguppi_xgpu_output_databuf_header/data()`
//! - don't size directly, use `sizeof_hpguppi_output_xgpu_block_t()`
typedef struct hpguppi_output_xgpu_block {
  char hdr[BLOCK_HDR_SIZE];
  char data[];
} hpguppi_output_xgpu_block_t;

typedef struct hpguppi_output_xgpu_databuf {
  hashpipe_databuf_t header;
  hashpipe_databuf_alignment padding; // Maintain data alignment
  hpguppi_output_xgpu_block_t block[N_XGPU_OUTPUT_BLOCKS]; //! this is dynamically sized, see comments
} hpguppi_output_xgpu_databuf_t;

static inline size_t hpguppi_output_xgpu_block_data_byte_size() {
  XGPUInfo xgpu_info = {0};
  xgpuInfo(&xgpu_info);
  return xgpu_info.matLength*sizeof(XGPU_OUTPUT_BLOCK_ELEMENT_T);
}

// #define sizeof(hpguppi_output_xgpu_block_t) (hpguppi_output_xgpu_block_data_byte_size()+BLOCK_HDR_SIZE)
static inline size_t sizeof_hpguppi_output_xgpu_block_t() {
    return hpguppi_output_xgpu_block_data_byte_size()+BLOCK_HDR_SIZE;
}

hashpipe_databuf_t *hpguppi_output_xgpu_databuf_create(int instance_id, int databuf_id);

#endif