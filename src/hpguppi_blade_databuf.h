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

#include "hpguppi_ata_blade_mode.h"

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

#endif // _HPGUPPI_BLADE_DATABUF_H
