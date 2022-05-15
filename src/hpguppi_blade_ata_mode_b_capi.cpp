#include <cassert>
#include <memory>

#include "blade/base.hh"
#include "blade/runner.hh"
#include "blade/pipelines/ata/mode_b.hh"

extern "C" {
#include "hpguppi_blade_ata_mode_b_capi.h"
}

using namespace Blade;
using namespace Blade::Pipelines::ATA;

using BladePipeline = ModeB<BLADE_ATA_MODE_B_OUTPUT_ELEMENT_T>;
static std::unique_ptr<Runner<BladePipeline>> runner;

bool blade_use_device(int device_id) {
    return SetCudaDevice(device_id) == Result::SUCCESS;
}

bool blade_ata_b_initialize(
    struct blade_ata_mode_b_config ata_b_config,
    size_t numberOfWorkers,
    struct blade_ata_mode_b_observation_meta* observationMeta,
    struct LonLatAlt* arrayReferencePosition,
    double* obs_phase_center_radecrad,
    double* beamCoordinates_radecrad,
    double* antennaPositions_xyz,
    double _Complex* antennaCalibrations
) {
    if (runner) {
        BL_FATAL("Can't initialize because Blade Runner is already initialized.");
        throw Result::ASSERTION_ERROR;
    }

    std::vector<XYZ> antennaPositions(ata_b_config.inputDims.NANTS);
    std::vector<RA_DEC> beamCoordinates(ata_b_config.beamformerBeams);
    std::vector<std::complex<double>> antennaCalibrationsCpp(
            ata_b_config.inputDims.NANTS*\
            ata_b_config.inputDims.NCHANS*\
            ata_b_config.inputDims.NPOLS);
    int i;
    for(i = 0; i < ata_b_config.inputDims.NANTS; i++){
        antennaPositions[i].X = antennaPositions_xyz[i*3 + 0];
        antennaPositions[i].Y = antennaPositions_xyz[i*3 + 1];
        antennaPositions[i].Z = antennaPositions_xyz[i*3 + 2];
    }
    for(i = 0; i < ata_b_config.beamformerBeams; i++){
        beamCoordinates[i].RA = beamCoordinates_radecrad[i*2 + 0];
        beamCoordinates[i].DEC = beamCoordinates_radecrad[i*2 + 1];
    }
    memcpy(antennaCalibrationsCpp.data(), antennaCalibrations,
            antennaCalibrationsCpp.size()*sizeof(antennaCalibrationsCpp[0]));

    BladePipeline::Config config = {
        .numberOfBeams = ata_b_config.inputDims.NBEAMS,
        .numberOfAntennas  = ata_b_config.inputDims.NANTS,
        .numberOfFrequencyChannels = ata_b_config.inputDims.NCHANS,
        .numberOfTimeSamples  = ata_b_config.inputDims.NTIME,
        .numberOfPolarizations  = ata_b_config.inputDims.NPOLS,

        .channelizerRate = ata_b_config.channelizerRate,
        .beamformerBeams = ata_b_config.beamformerBeams,

        .rfFrequencyHz = observationMeta->rfFrequencyHz,
        .channelBandwidthHz = observationMeta->channelBandwidthHz,
        .totalBandwidthHz = observationMeta->totalBandwidthHz,
        .frequencyStartIndex = observationMeta->frequencyStartIndex,
        .referenceAntennaIndex = observationMeta->referenceAntennaIndex,
        .arrayReferencePosition = {
            .LON = arrayReferencePosition->LON,
            .LAT = arrayReferencePosition->LAT,
            .ALT = arrayReferencePosition->ALT
        },
        .boresightCoordinate = {
            .RA = obs_phase_center_radecrad[0],
            .DEC = obs_phase_center_radecrad[1]
        },
        .antennaPositions = antennaPositions,
        .antennaCalibrations = antennaCalibrationsCpp,
        .beamCoordinates = beamCoordinates,

        .outputMemWidth = ata_b_config.outputMemWidth,
        .outputMemPad = ata_b_config.outputMemPad,

        .castBlockSize = ata_b_config.castBlockSize,
        .channelizerBlockSize = ata_b_config.channelizerBlockSize,
        .beamformerBlockSize = ata_b_config.beamformerBlockSize,
    };

    //config.antennaCalibrations.resize(
    //    config.numberOfAntennas *
    //    config.numberOfFrequencyChannels *
    //    config.channelizerRate *
    //    config.numberOfPolarizations
    //);
    
    runner = Runner<BladePipeline>::New(numberOfWorkers, config);

    return true;
}

void blade_ata_b_terminate() {
    if (!runner) {
        BL_FATAL("Can't terminate because Blade Runner isn't initialized.");
        throw Result::ASSERTION_ERROR;
    }
    runner.reset();
}

size_t blade_ata_b_get_input_size() {
    assert(runner);
    return runner->getWorker().getInputSize();
}

size_t blade_ata_b_get_output_size() {
    assert(runner);
    return runner->getWorker().getOutputSize();
}

bool blade_pin_memory(void* buffer, size_t size) {
    return Memory::PageLock(Vector<Device::CPU, I8>(buffer, size)) == Result::SUCCESS;
}

bool blade_ata_b_enqueue(void* input_ptr, void* output_ptr, size_t id, double time_mjd, double dut1) {
    assert(runner);
    return runner->enqueue([&](auto& worker){
        auto input = Vector<Device::CPU, CI8>(input_ptr, worker.getInputSize());
        auto output = Vector<Device::CPU, BLADE_ATA_MODE_B_OUTPUT_ELEMENT_T>(output_ptr, worker.getOutputSize());

        worker.run(time_mjd, dut1, input, output);

        return id;
    });
}

bool blade_ata_b_dequeue(size_t* id) {
    assert(runner);
    return runner->dequeue(id);
}
