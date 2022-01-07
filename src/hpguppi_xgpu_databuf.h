/* hpguppi_databuf.h
 *
 * Defines shared mem structure for data passing.
 */
#ifndef _HPGUPPI_XGPU_DATABUF_H
#define _HPGUPPI_XGPU_DATABUF_H
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

static inline hpguppi_input_xgpu_databuf_t *hpguppi_input_xgpu_databuf_attach(int instance_id, int databuf_id)
{
    return (hpguppi_input_xgpu_databuf_t *)hashpipe_databuf_attach(instance_id, databuf_id);
}

static inline int hpguppi_input_xgpu_databuf_detach(hpguppi_input_xgpu_databuf_t *d)
{
    return hashpipe_databuf_detach((hashpipe_databuf_t *)d);
}

static inline void hpguppi_input_xgpu_databuf_clear(hpguppi_input_xgpu_databuf_t *d)
{
    hashpipe_databuf_clear((hashpipe_databuf_t *)d);
}

static inline int hpguppi_input_xgpu_databuf_block_status(hpguppi_input_xgpu_databuf_t *d, int block_id)
{
    return hashpipe_databuf_block_status((hashpipe_databuf_t *)d, block_id);
}

static inline int hpguppi_input_xgpu_databuf_total_status(hpguppi_input_xgpu_databuf_t *d)
{
    return hashpipe_databuf_total_status((hashpipe_databuf_t *)d);
}

static inline int hpguppi_input_xgpu_databuf_wait_free_timeout(
    hpguppi_input_xgpu_databuf_t *d, int block_id, struct timespec *timeout)
{
    return hashpipe_databuf_wait_free_timeout((hashpipe_databuf_t *)d,
        block_id, timeout);
}

static inline int hpguppi_input_xgpu_databuf_wait_free(hpguppi_input_xgpu_databuf_t *d, int block_id)
{
    return hashpipe_databuf_wait_free((hashpipe_databuf_t *)d, block_id);
}

static inline int hpguppi_input_xgpu_databuf_busywait_free(hpguppi_input_xgpu_databuf_t *d, int block_id)
{
    return hashpipe_databuf_busywait_free((hashpipe_databuf_t *)d, block_id);
}

static inline int hpguppi_input_xgpu_databuf_wait_filled_timeout(
    hpguppi_input_xgpu_databuf_t *d, int block_id, struct timespec *timeout)
{
    return hashpipe_databuf_wait_filled_timeout((hashpipe_databuf_t *)d,
        block_id, timeout);
}

static inline int hpguppi_input_xgpu_databuf_wait_filled(hpguppi_input_xgpu_databuf_t *d, int block_id)
{
    return hashpipe_databuf_wait_filled((hashpipe_databuf_t *)d, block_id);
}

static inline int hpguppi_input_xgpu_databuf_busywait_filled(hpguppi_input_xgpu_databuf_t *d, int block_id)
{
    return hashpipe_databuf_busywait_filled((hashpipe_databuf_t *)d, block_id);
}

static inline int hpguppi_input_xgpu_databuf_set_free(hpguppi_input_xgpu_databuf_t *d, int block_id)
{
    return hashpipe_databuf_set_free((hashpipe_databuf_t *)d, block_id);
}

static inline int hpguppi_input_xgpu_databuf_set_filled(hpguppi_input_xgpu_databuf_t *d, int block_id)
{
    return hashpipe_databuf_set_filled((hashpipe_databuf_t *)d, block_id);
}

static inline char *hpguppi_xgpu_databuf_header(struct hpguppi_input_xgpu_databuf *d, int block_id) {
    if(block_id < 0 || d->header.n_block < block_id) {
        hashpipe_error(__FUNCTION__,
            "block_id %s out of range [0, %d)",
            block_id, d->header.n_block);
        return NULL;
    } else {
        return d->block[block_id].hdr;
    }
}

static inline char *hpguppi_xgpu_databuf_data(struct hpguppi_input_xgpu_databuf *d, int block_id) {
    if(block_id < 0 || d->header.n_block < block_id) {
        hashpipe_error(__FUNCTION__,
            "block_id %s out of range [0, %d)",
            block_id, d->header.n_block);
        return NULL;
    } else {
        return d->block[block_id].data;
    }
}

#endif