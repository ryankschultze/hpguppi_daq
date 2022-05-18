#ifndef BLADE_ATA_STRUCTS_H
#define BLADE_ATA_STRUCTS_H

struct blade_ata_input_dims{
    uint32_t NBEAMS;
    uint32_t NANTS;
    uint32_t NCHANS;
    uint32_t NTIME;
    uint32_t NPOLS;
};

struct blade_ata_observation_meta{
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

#endif //BLADE_ATA_STRUCTS_H