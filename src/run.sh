#!/bin/bash
# Set high performance mode
for i in `ls /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor`; do echo performance > $i; done

perf=
if [ "$1" = 'perf' ]
then
  perf=perf
  shift
fi

# Set mtu
ifconfig enp134s0d1 mtu 9000

# Kernel buffer sizes
sysctl net.core.rmem_max=1073741824 # 1G 134217728 #128MB
sysctl net.core.rmem_default=1073741824 # 1G 134217728 #128MB
# sysctl net.core.rmem_max=65536 # 1G 134217728 #128MB
# sysctl net.core.rmem_default=65536 # 1G 134217728 #128MB

# Kill packets before the IP stack, applicable to hashpipe_pktsock
iptables -t raw -A PREROUTING -i enp134s0d1 -p udp -j DROP

# Set interrupt coalescing
ethtool -C enp134s0d1 adaptive-rx on
ethtool -C enp134s0d1 rx-frames 8
ethtool -C enp134s0d1 rx-usecs 0

# Set ring sizes to max
ethtool -G enp134s0d1 rx 8192

numactl --cpunodebind=1 --membind=1 \
  $perf \
hashpipe -p /usr/local/lib/hpguppi_daq.so -I 0 \
-o BINDHOST=enp134s0d1 \
-o BINDPORT=10000 \
-c 9 hpguppi_atasnap_pktsock_thread -m 57344 hpguppi_atasnap_pkt_to_FTP_transpose -c 11 hpguppi_atasnap_rawdisk_thread #null_output_thread
# -c 9 hpguppi_atasnap_pktsock_thread -c 10 hpguppi_rawdisk_only_thread #null_output_thread

# echo "Starting hashpipe REDIS Gateway"
# # hashpipe_redis_gateway.rb -s ${REDISHOST:-redishost}# -D hashpipe -g `hostname -s` -i 0 -f &
# hashpipe_redis_gateway.rb -s ${REDISHOST:-redishost} -f & # -D hashpipe -g `hostname -s` -i 0 -f &
