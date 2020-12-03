#!/bin/bash
export PSRHOME="/home/obsuser/"

# quietsnaps="frb-snap8-pi" # frb-snap7-pi frb-snap9-pi"
# quietsnaps="frb-snap9-pi"
snaps="frb-snap8-pi"
# snapsA="frb-snap5-pi"
# snapsB="frb-snap7-pi"
 
# Stop ethernet
# /home/obsuser/nov_observing/stop_snaps_eth.py $quietsnaps

# Start ethernet
# /home/obsuser/nov_observing/start_snaps_eth.py $snaps

# /home/obsuser/nov_observing/sync_snaps.py $snaps #resets the pkt_idx of each snap

/opt/mnt/bin/hashpipe_check_status -k 	OBSSTART  -i 	0
/opt/mnt/bin/hashpipe_check_status -k 	OBSSTOP  	-i 	0

echo "Setup the Hashpipe Status keys"
/home/obsuser/nov_observing/populate_meta.py /home/sonata/src/hpguppi_daq/ataconfig.yml -s $snaps -i
# /home/obsuser/nov_observing/populate_meta.py ~/mydonsol/hpguppi_daq/src/ataconfig-port10000.yml -s $snapsA -i
# /home/obsuser/nov_observing/populate_meta.py ~/mydonsol/hpguppi_daq/src/ataconfig-port9000.yml -s $snapsB -i

echo "Start recording observation with /home/obsuser/nov_observing/start_record_in_x.py"