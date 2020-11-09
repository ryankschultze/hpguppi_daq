#!/bin/bash
# Set high performance mode
for i in `ls /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor`; do echo performance > $i; done

# Set mtu
ifconfig enp134s0d1 mtu 9000

# Kernel buffer sizes
sysctl net.core.rmem_max=1073741824 #1G
sysctl net.core.rmem_default=1073741824 #1G

# Kill packets before the IP stack, applicable to hashpipe_pktsock
iptables -t raw -A PREROUTING -i enp134s0d1 -p udp -j DROP

# Set interrupt coalescing
ethtool -C enp134s0d1 adaptive-rx on
ethtool -C enp134s0d1 rx-frames 8
ethtool -C enp134s0d1 rx-usecs 0

# Set ring sizes to max
ethtool -G enp134s0d1 rx 8192

hashpipe -p /usr/local/lib/hpguppi_daq.so -I 0 \
-o BINDHOST=enp134s0d1 \
-o BINDPORT=4015 \
-o DATADIR=/mnt/buf0 \
-o DESTIP=10.11.1.156 \
-o DESTDIR=/mnt/buf0 \
-o FENCHAN=4096 \
-o NANTS=1 \
-o NPOL=2 \
-o NBITS=4  \
-o NSTRM=1 \
-o PKTNTIME=16 \
-o OBSSTART=0 \
-o OBSSTOP=67108864 \
-o SCHAN=0 \
-o CHAN_BW=0.25 \
-o OBSBW=128 \
-o PROJID=dmp \
-o PKTFMT=ATASNAPV \
-o OBSFREQ=1420 \
-o OBSNCHAN=256 -c 9 hpguppi_atasnap_pktsock_thread -c 10 null_output_thread #hpguppi_rawdisk_only_thread