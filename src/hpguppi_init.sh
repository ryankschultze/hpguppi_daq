#!/bin/bash

# Add directory containing this script to PATH
PATH="$(dirname $0):${PATH}"

PREFIX_EXEC="/opt/mnt/bin/"
PREFIX_LIB="/opt/mnt/lib/"

hostname=`hostname -s`

cores_per_cpu=`grep cpu.cores /proc/cpuinfo | awk 'NR==1{print $NF}'`

case $cores_per_cpu in
  6)
    # We use a dual hexa-core CPU system with everything on the second socket (CPUs 6-11),
    # with CPU 8 being reserved for NIC IRQs.  This leads to the following CPU mask:
    #
    # 1111 1100 0000 0000
    # 5432 1098 7654 3210
    # ---- ---- ---- ----
    # 0000 1110 1100 0000 = 0x0ec0
    NUMABIND0=1
    DATADIR0=/mnt/buf0
    NUMABIND1=1
    DATADIR0=/mnt/buf0
    NET0CPU=9
    OUT0CPU=10
    TSP0MSK=63
    BINDHOST0=eth4
    NET1CPU=6
    OUT1CPU=7
    TSP1MSK=63
    BINDHOST1=eth5
    ;;
  8)
    # We use a dual octa-core CPU system with everything on the second socket (CPUs 8-15),
    # with CPU 8 being reserved for NIC IRQs.  This leads to the following CPU mask:
    #
    # 1111 1100 0000 0000
    # 5432 1098 7654 3210
    # ---- ---- ---- ----
    # 1111 1110 0000 0000 = 0xfe00
    NUMABIND0=1
    DATADIR0=/mnt/buf0
    NUMABIND1=1
    DATADIR0=/mnt/buf0
    NET0CPU=8
    OUT0CPU=9
    TPS0MSK=30720 #59392 #30720 #61440 #57344 #49152
    BINDHOST0=ens6d1
    NET1CPU=13
    OUT1CPU=14
    TPS1MSK=61440
    BINDHOST1=ens6d0
    ;;
  16)
    # We use a dual 16-core CPU system with each socket having dedicated
    # NICs and mounts.  This leads to the following CPU mask
    # (0s for instance 0, 1s for instance 1):
    #
    # 3322 2222 2222 2222 1111 1100 0000 0000
    # 1098 7654 3210 9876 5432 1098 7654 3210
    # ---- ---- ---- ---- ---- ---- ---- ----
    # 1111 1111 1111 1111 0000 0000 0000 0000 = 0xffff0000
    NUMABIND0=1
    DATADIR0=/mnt/buf0
    NET0CPU=17
    OUT0CPU=18
    TPS0MSK=133693440
    BINDHOST0=enp97s0f1
    NUMABIND1=0
    DATADIR1=/mnt/buf1
    NET1CPU=1
    OUT1CPU=2
    TPS1MSK=2040
    BINDHOST1=enp225s0f1
    ;;
  *)
    echo "$cores_per_cpu cores per cpu is not yet supported by $0"
    exit 1
    ;;
esac

instances=(
  # Setup parameters for two instances.
  #
  #
  #              cpu                bind          NET       OUT     log
  #  dir        socket    mask      host          CPU       CPU     dir
  "$DATADIR0  $NUMABIND0 $TPS0MSK  $BINDHOST0  $NET0CPU  $OUT0CPU  /home/sonata/logs"  # Instance 0
  "$DATADIR1  $NUMABIND1 $TPS1MSK  $BINDHOST1  $NET1CPU  $OUT1CPU  /home/sonata/logs"  # Instance 1
)

function init() {
  instance=$1
  dir=$2
  numabind=$3
  mask=$4
  bindhost=$5
  netcpu=$6
  outcpu=$7
  workdir=$8
  shift 8

  if [ $numabind -lt 0 ]
  then
    numabind=$instance
  fi

  if [ -z "${dir}" ]
  then
    echo "Invalid instance number '${instance:-[unspecified]}' (ignored)"
    return 1
  #elif [ "${dir}" == "/buf0" ]
  # elif [ "${dir%/buf?}" != "${dir}" ]
  # then
  #   # Don't want to output messages to NVMe directory
  #   echo "setting workdir to /tmp"
  #   workdir=/tmp
  fi

  if [ -z "$outcpu" ]
  then
    echo "Invalid configuration for host ${hostname} instance ${instance} (ignored)"
    return 1
  fi

  if ! cd "${workdir}"
  then
    echo "Invalid working directory ${workdir} (ignored)"
    return 1
  fi

  if test -f "${workdir}/${hostname}.$instance.out"; then
    echo "Trimming log ${workdir}/${hostname}.$instance.out"
    tail -n 100000 ${workdir}/${hostname}.$instance.out > tmp.out
    mv tmp.out ${workdir}/${hostname}.$instance.out
  fi
  if test -f "${workdir}/${hostname}.$instance.err"; then
    echo "Trimming log ${workdir}/${hostname}.$instance.err"
    tail -n 100000 ${workdir}/${hostname}.$instance.err > tmp.err
    mv tmp.err ${workdir}/${hostname}.$instance.err
  fi

  echo numactl --cpunodebind=$numabind --membind=$numabind \
  $perf \
  ${PREFIX_EXEC}hashpipe -p ${hpguppi_plugin:-hpguppi_daq} -I $instance \
    -o BINDHOST=${bindhost}${vlan} \
    ${bindport:+-o BINDPORT=$bindport} \
    -o DATADIR=$dir \
    ${@} \
    -c $netcpu $net_thread \
    -m $mask $work_thread \
    -c $outcpu $out_thread

  numactl --cpunodebind=$numabind --membind=$numabind \
  $perf \
  ${PREFIX_EXEC}hashpipe -p ${hpguppi_plugin:-hpguppi_daq} -I $instance \
    -o BINDHOST=${bindhost}${vlan} \
    ${bindport:+-o BINDPORT=$bindport} \
    -o DATADIR=$dir \
    ${@} \
    -c $netcpu $net_thread \
    -m $mask $work_thread \
    -c $outcpu $out_thread < /dev/null 1>> "${workdir}/${hostname}.$instance.out" 2>> "${workdir}/${hostname}.$instance.err" &
}

perf=
redis_sync_key=sync_time
net_thread=hpguppi_net_thread
out_thread=hpguppi_rawdisk_thread
work_thread=
bindport=60000
vlan=
use_fifo=yes
options=
instance_ports=(
  ""
  ""
)

if [ "$1" = 'perf' ]
then
  perf="perf record"
  shift
fi

if [ "$1" = 'fakefake' ]
then
  net_thread=hpguppi_fake_net_thread
  out_thread=null_output_thread
  shift
elif [ "$1" = 'fake' ]
then
  net_thread=hpguppi_fake_net_thread
  shift
elif [ "$1" = 'mb1' ]
then
  redis_sync_key=s6_mcount_0
  net_thread=hpguppi_mb1_net_thread
  bindport=$((0x5336)) # ASCII 0x53 0x36 is "S6" (for SEREDIP6)
  shift
elif [ "$1" = 'mb128ch' ]
then
  redis_sync_key=s6_mcount_0
  net_thread=hpguppi_mb128ch_net_thread
  bindport=$((0x4D42)) # ASCII 0x4D 0x42 is "MB" (for MultiBeam)
  shift
elif [ "$1" = 'pksuwl' ]
then
  redis_sync_key=
  use_fifo=no
  hpguppi_plugin=/home/davidm/local/src/hpguppi_daq/src/.libs/hpguppi_daq.so
  #net_thread=hpguppi_pksuwl_net_thread
  net_thread="hpguppi_ibvpkt_thread -c 11 hpguppi_pksuwl_vdif_thread"
  out_thread=hpguppi_fildisk_only_thread
  # BINDPORT must be passed as a separate "-o BINDPORT=1234" command line
  # option by caller (e.g. hpguppi_pksuwl_init.sh)
  bindport=
  vlan='.2'
  #options="-o IBVPKTSZ=42,32,8198" # Not 802.1Q tagged VLAN
  options="-o IBVPKTSZ=46,32,8198" # 802.1Q tagged VLAN
  # Handle special case of specifying disk thread as last command line param
  case "${@:$#}" in
    *_thread)
      # Use last param as out_thread
      out_thread="${@:$#}"
      # Remove last param from $@
      set -- "${@:1:$(($#-1))}"
      ;;
  esac
  shift
elif [ "$1" = 'meerkat' ]
then
  redis_sync_key=
  use_fifo=no
  #net_thread=hpguppi_meerkat_net_thread
  #net_thread=hpguppi_meerkat_pkt_thread
  net_thread="hpguppi_ibvpkt_thread -c 11 hpguppi_meerkat_spead_thread"
  options="-o IBVPKTSZ=42,96,1024"
  bindport=7148
  instances[0]="${instances[0]/datax/buf0}"
  instances[1]="${instances[1]/datax2/buf1}"
  # For initial testing...
  out_thread=null_output_thread
  #out_thread=hpguppi_rawdisk_only_thread
  shift
elif [ "$1" = 'atasnap' ]
then
  use_fifo=no
  # net_thread="hpguppi_ibvpkt_thread -c 11 hpguppi_atasnap_voltage_thread"
  # net_thread="hpguppi_net_thread -c 11 hpguppi_atasnap_voltage_thread"
  # net_thread="hpguppi_net_thread"

  net_thread=hpguppi_atasnap_pktsock_thread
  work_thread=hpguppi_atasnap_pkt_to_FTP_transpose
  # options="-o IBVPKTSZ=42,8,8192"
  instances[0]="${instances[0]/eth4/ens6d1}"
  instances[1]="${instances[1]/eth5/ens6d1}"
  bindport=0
  instance_ports[0]=10000
  instance_ports[1]=10000
  # instances[1]="${instances[0]/datax/mnt/buf0}"
  # instances[1]="${instances[1]/datax2/mnt/buf1}"
  # For initial testing...
  # hpguppi_plugin=/homelocal/sonata/davidm/src/hpguppi_daq/src/.libs/hpguppi_daq.so
  hpguppi_plugin=${PREFIX_LIB}hpguppi_daq.so
  # out_thread=null_output_thread
  # out_thread=hpguppi_rawdisk_only_thread
  out_thread=hpguppi_atasnap_rawdisk_thread
  shift
elif echo "$1" | grep -q 'thread'
then
  out_thread="$1"
  shift
fi

# Exit if no instance id is given
if [ -z "$1" ]
then
  echo "Usage: $(basename $0) [mb1|mb128ch|pksuwl|meerkat|fake|fakefake] INSTANCE_ID [...] [OPTIONS]"
  exit 1
fi

echo using net_thread $net_thread
echo using out_thread $out_thread

instance_ids=''

while [ -n "$1" ]
do
  # Break on first option character
  case $1 in
    -*) break;;
  esac

  instance_ids="${instance_ids} ${1}"
  shift
done

# Get sync time from redis
if [ -n "${redis_sync_key}" ]
then
  sync_time=`redis-cli -h redishost get $redis_sync_key | tr -d '"'`
else
  sync_time=0
fi

for instidx in $instance_ids
do
  fifo="/tmp/hpguppi_daq_control/$instidx"
  args="${instances[$instidx]}"
  if [ -n "${args}" ]
  then
    echo Starting instance $hostname/$instidx
    bindport=${instance_ports[$instidx]}
    if init $instidx $args $options "${@}"
    then
      echo Instance $hostname/$instidx pid $!
      ${PREFIX_EXEC}hashpipe_check_status -I $instidx -k SYNCTIME -i ${sync_time:-0}

      test "$use_fifo" == "yes" || continue

      # Loop until control fifo exists, time out after 3 seconds
      timeout=3
      while [ $timeout -gt 0 -a ! -e $fifo ]
      do
        timeout=$((timeout-1))
        sleep 1
      done

      # If FIFO exists, use it to "start" the instance
      # (i.e. put it in the armed state)
      if [ -e $fifo ]
      then
        echo start > $fifo
      else
        echo "Error: control FIFO not created within 3 seconds"
      fi
    fi
  else
    echo Instance $instidx not defined for host $hostname
  fi
done
