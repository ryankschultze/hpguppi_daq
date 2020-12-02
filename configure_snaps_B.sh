#!/bin/bash

export PSRHOME="/home/obsuser"
fpgfile="/home/obsuser/src/ata_snap/snap_adc5g_feng_rpi/outputs/snap_adc5g_feng_rpi_2020-11-26_1432.fpg"

# configuration file to use
cfg="/homelocal/sonata/mydonsol/hpguppi_daq/ataconfig-port9000.yml"
cmd=/home/obsuser/miniconda3/envs/ATAobs/bin/snap_feng_init.py

# snap physical ids, eg 5 -> frb-snap5-pi
snapids=( 7 )

echo "-------"
echo "0) Resetting the OBSSTART/OBSSTOP before reconfiguring"
echo /home/obsuser/nov_observing/start_record_in_x.py -r
/home/obsuser/nov_observing/start_record_in_x.py -r

echo "-------"
echo "1) Reprogramming/configuring snaps"
#for snapid in {1,2,3,4,5,6,7,8,9,10,11,12}
snaps=""
ifeng=0
for snapid in "${snapids[@]}"
do
    echo $cmd frb-snap$snapid-pi $fpgfile $cfg -i $ifeng -s --eth_volt --skipprog
    $cmd frb-snap$snapid-pi $fpgfile $cfg -i $ifeng -s --eth_volt --skipprog
    snaps="${snaps}frb-snap${snapid}-pi "
    ((ifeng=ifeng+1))
done

echo "-------"
echo "2) Syncing snaps"
echo /home/obsuser/nov_observing/sync_snaps.py ${snaps}
/home/obsuser/nov_observing/sync_snaps.py ${snaps}

echo
echo "-------"
echo "3) Populating REDIS meta data"
echo /homelocal/sonata/mydonsol/hpguppi_daq/populate_meta.py ${cfg} -s ${snaps} -I -i 1
/homelocal/sonata/mydonsol/hpguppi_daq/populate_meta.py ${cfg} -s ${snaps} -I -i 1
