/* hpguppi_databuf.h
 *
 * Defines shared mem structure for data passing.
 */
#ifndef _HPGUPPI_BLADE_DATABUF_H
#define _HPGUPPI_BLADE_DATABUF_H
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

#include "blade/pipelines/ata/mode_b_config.h"

#define BLADE_BLOCK_DATA_SIZE (BLADE_OUTPUT_NBEAM *\
                               BLADE_ANT_NCHAN *\
                               BLADE_NTIME *\
                               BLADE_NPOL *\
                               BLADE_OUTPUT_NCOMPLEX_BYTES)

typedef struct {
  char hdr[BLOCK_HDR_SIZE];
  char data[BLADE_BLOCK_DATA_SIZE];
} hpguppi_blade_block_t;

// Used to pad after hashpipe_databuf_t to maintain data alignment
typedef uint8_t hashpipe_databuf_alignment[
  ALIGNMENT_SIZE - (sizeof(hashpipe_databuf_t)%ALIGNMENT_SIZE)
];

typedef struct  {
  hashpipe_databuf_t header;
  hashpipe_databuf_alignment padding; // Maintain data alignment
  hpguppi_blade_block_t block[N_INPUT_BLOCKS];
} hpguppi_blade_output_databuf_t;

/*
 * BLADE_OUTPUT BUFFER FUNCTIONS
 */

hashpipe_databuf_t *hpguppi_blade_output_databuf_create(int instance_id, int databuf_id);

static inline hpguppi_blade_output_databuf_t *hpguppi_blade_output_databuf_attach(int instance_id, int databuf_id)
{
    return (hpguppi_blade_output_databuf_t *)hashpipe_databuf_attach(instance_id, databuf_id);
}

static inline int hpguppi_blade_output_databuf_detach(hpguppi_blade_output_databuf_t *d)
{
    return hashpipe_databuf_detach((hashpipe_databuf_t *)d);
}

static inline void hpguppi_blade_output_databuf_clear(hpguppi_blade_output_databuf_t *d)
{
    hashpipe_databuf_clear((hashpipe_databuf_t *)d);
}

static inline int hpguppi_blade_output_databuf_block_status(hpguppi_blade_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_block_status((hashpipe_databuf_t *)d, block_id);
}

static inline int hpguppi_blade_output_databuf_total_status(hpguppi_blade_output_databuf_t *d)
{
    return hashpipe_databuf_total_status((hashpipe_databuf_t *)d);
}

static inline int hpguppi_blade_output_databuf_wait_free_timeout(
    hpguppi_blade_output_databuf_t *d, int block_id, struct timespec *timeout)
{
    return hashpipe_databuf_wait_free_timeout((hashpipe_databuf_t *)d,
        block_id, timeout);
}

static inline int hpguppi_blade_output_databuf_wait_free(hpguppi_blade_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_wait_free((hashpipe_databuf_t *)d, block_id);
}

static inline int hpguppi_blade_output_databuf_busywait_free(hpguppi_blade_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_busywait_free((hashpipe_databuf_t *)d, block_id);
}

static inline int hpguppi_blade_output_databuf_wait_filled_timeout(
    hpguppi_blade_output_databuf_t *d, int block_id, struct timespec *timeout)
{
    return hashpipe_databuf_wait_filled_timeout((hashpipe_databuf_t *)d,
        block_id, timeout);
}

static inline int hpguppi_blade_output_databuf_wait_filled(hpguppi_blade_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_wait_filled((hashpipe_databuf_t *)d, block_id);
}

static inline int hpguppi_blade_output_databuf_busywait_filled(hpguppi_blade_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_busywait_filled((hashpipe_databuf_t *)d, block_id);
}

static inline int hpguppi_blade_output_databuf_set_free(hpguppi_blade_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_set_free((hashpipe_databuf_t *)d, block_id);
}

static inline int hpguppi_blade_output_databuf_set_filled(hpguppi_blade_output_databuf_t *d, int block_id)
{
    return hashpipe_databuf_set_filled((hashpipe_databuf_t *)d, block_id);
}

static inline char *hpguppi_blade_databuf_header(hpguppi_blade_output_databuf_t *d, int block_id) {
    if(block_id < 0 || d->header.n_block < block_id) {
        hashpipe_error(__FUNCTION__,
            "block_id %s out of range [0, %d)",
            block_id, d->header.n_block);
        return NULL;
    } else {
        return d->block[block_id].hdr;
    }
}

static inline char *hpguppi_blade_databuf_data(hpguppi_blade_output_databuf_t *d, int block_id) {
    if(block_id < 0 || d->header.n_block < block_id) {
        hashpipe_error(__FUNCTION__,
            "block_id %s out of range [0, %d)",
            block_id, d->header.n_block);
        return NULL;
    } else {
        return d->block[block_id].data;
    }
}

#endif // _HPGUPPI_BLADE_DATABUF_H
