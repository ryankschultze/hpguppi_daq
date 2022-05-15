#ifndef MODEB_H
#define MODEB_H

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <complex.h>

#include "hpguppi_blade_ata_mode_b_config.h"

struct blade_ata_mode_b_input_dims{
    uint32_t NBEAMS;
    uint32_t NANTS;
    uint32_t NCHANS;
    uint32_t NTIME;
    uint32_t NPOLS;
};

struct blade_ata_mode_b_observation_meta{
    double rfFrequencyHz; // OBSFREQ
    double channelBandwidthHz; // CHAN_BW
    double totalBandwidthHz; // CHAN_BW*FENCHAN
    uint64_t frequencyStartIndex; // SCHAN
    uint64_t referenceAntennaIndex; // ???
};

struct LonLatAlt{
    double LON;
    double LAT;
    double ALT;
};

struct blade_ata_mode_b_config {
    struct blade_ata_mode_b_input_dims inputDims;
    uint32_t channelizerRate;
    uint32_t beamformerBeams;

    uint32_t outputMemWidth;
    uint32_t outputMemPad;

    uint32_t castBlockSize;
    uint32_t channelizerBlockSize;
    uint32_t beamformerBlockSize;
};

const struct blade_ata_mode_b_config BLADE_ATA_MODE_B_CONFIG = {
    { // .inputDims
        1, // .NBEAMS
        BLADE_ATA_MODE_B_INPUT_NANT, // .NANTS 
        BLADE_ATA_MODE_B_ANT_NCHAN, // .NCHANS
        BLADE_ATA_MODE_B_NTIME, // .NTIME 
        BLADE_ATA_MODE_B_NPOL, // .NPOLS 
    }, // .inputDims
    BLADE_ATA_MODE_B_CHANNELIZER_RATE, // .channelizerRate
    BLADE_ATA_MODE_B_OUTPUT_NBEAM, // .beamformerBeams

    BLADE_ATA_MODE_B_OUTPUT_MEMCPY2D_WIDTH, // .outputMemWidth
    BLADE_ATA_MODE_B_OUTPUT_MEMCPY2D_PAD, // .outputMemPad

    512, // .castBlockSize
    512, // .channelizerBlockSize
    512  // .beamformerBlockSize
};

bool blade_use_device(int device_id);
bool blade_ata_b_initialize(
    struct blade_ata_mode_b_config ata_b_config,
    size_t numberOfWorkers,
    struct blade_ata_mode_b_observation_meta* observationMeta,
    struct LonLatAlt* arrayReferencePosition,
    double* obs_phase_center_radecrad,
    double* beamCoordinates_radecrad,
    double* antennaPositions_xyz,
    double _Complex* antennaCalibrations
);
size_t blade_ata_b_get_input_size();
size_t blade_ata_b_get_output_size();
bool blade_pin_memory(void* buffer, size_t size);
bool blade_ata_b_set_antenna_positions(void* xyz_positions, bool block);
bool blade_ata_b_set_antenna_calibrations(void* calibrations, bool block);
bool blade_ata_b_set_beam_coordinates(void* coordinates, bool block);
bool blade_ata_b_set_boresight_coordinates(void* coordinate, bool block);
bool blade_ata_b_enqueue(void* input_ptr, void* output_ptr, size_t id, double time_mjd, double dut1);
bool blade_ata_b_dequeue(size_t* id);
void blade_ata_b_terminate();

#endif
