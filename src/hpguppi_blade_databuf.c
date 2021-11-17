/* hpguppi_blade_databuf.c
 *
 * Routines for creating and accessing BLADE_output_data transfer
 * buffer in shared memory.
 */

#include "hpguppi_blade_databuf.h"

hashpipe_databuf_t *hpguppi_blade_output_databuf_create(int instance_id, int databuf_id)
{
    int i;

    /* Calc databuf sizes */
    size_t header_size = sizeof(hashpipe_databuf_t)
                       + sizeof(hashpipe_databuf_alignment);
    size_t block_size  = sizeof(hpguppi_blade_block_t);
    int    n_block = N_INPUT_BLOCKS;

    hpguppi_blade_output_databuf_t * d = (hpguppi_blade_output_databuf_t *)
        hashpipe_databuf_create(
            instance_id, databuf_id, header_size, block_size, n_block);

    if(!d) {
      return NULL;
    }

    /* Zero out blocks */
    for(i=0; i<n_block; i++) {
      memset(&(d->block[i]), 0, sizeof(hpguppi_blade_block_t));
    }

    /* Init headers of each databuf block */
    char end_key[81];
    memset(end_key, ' ', 80);
    strncpy(end_key, "END", 3);
    end_key[80]='\0';
    for (i=0; i<n_block; i++) { 
        memcpy(hpguppi_blade_databuf_header(d,i), end_key, 80); 
    }

    return (hashpipe_databuf_t *)d;
}