#!/home/obsuser/miniconda3/envs/ATAobs/bin/python
import redis
import yaml
import argparse
import socket
import numpy as np
import sys
from string import Template

from SNAPobs import snap_defaults, snap_config
from ATATools import ata_control

NBITS              = 4
NPOL               = 2
REDISHOST          = 'redishost'
N_TIMES_PER_PKT    = 16
MAX_CHANS_PER_PKT  = 8*8192 // (2*NBITS) // N_TIMES_PER_PKT // 2
FENCHAN            = snap_defaults.nchan
BW                 = snap_defaults.bw * 1e6 # in Hz
FOFF               = snap_defaults.bw / snap_defaults.nchan # in Mhz
TBIN               = 1/(BW/FENCHAN) # in seconds

PROJID             = 'dmpauto_wael'
DATADIR            = '/mnt/buf0'
BACKEND            = 'GUPPI'
BANK               = '.'
PKTFMT             = 'ATASNAPV'
SNAPPAT            = 'frb-snap,pi'
BINDHOST           = 'enp134s0d1'

REDISSETGW = Template('hashpipe://${host}/${inst}/set')
REDISSET = 'hashpipe:///set'


ATA_SNAP_TAB = snap_config.get_ata_snap_tab()


def get_snap_mapping(snap_hosts, ignore_control=False):
    """
    1) Get the centre frequency for each of the snap_hosts
    from the REST/control machine, given what tunings are fed to 
    what snaps
    2) Get the antenna name associated with the snap

    An overly complicated function for what it does really
    """
    if type(snap_hosts) != list:
        raise RuntimeError("Please input a list")
    if not all(snap in list(ATA_SNAP_TAB.snap_hostname) for snap in snap_hosts):
        raise RuntimeError("Not all snaps (%s) are provided in the config table (%s)",
                snap_hosts, ATA_SNAP_TAB.snap_hostname)

    obs_ant_tab = ATA_SNAP_TAB[ATA_SNAP_TAB.snap_hostname.isin(snap_hosts)]
    los = np.unique(obs_ant_tab.LO)

    retdict_skyfreq = {}
    for lo in los:
        hostname_sub_list = list(obs_ant_tab[obs_ant_tab.LO == lo].snap_hostname)
        if ignore_control:
            skyfreq = 1400
        else:
            skyfreq = ata_control.get_sky_freq(lo=lo)
        tmp_dict = {snap:freq for snap,freq in zip(hostname_sub_list, 
            [skyfreq]*len(hostname_sub_list))}
        retdict_skyfreq.update(tmp_dict)

    retdict_antname = {i.snap_hostname:i.ANT_name for i in obs_ant_tab.itertuples()}
    return retdict_skyfreq, retdict_antname


def get_channel_selection(dests, start_chan, n_chans_per_dest):
    mapping = {}
    for dn,d in enumerate(dests):
        mapping[d] = list(range(start_chan + dn*n_chans_per_dest,
            start_chan + (dn+1)*n_chans_per_dest))
    return mapping


def main():
    parser = argparse.ArgumentParser(description='Program to populate'
            'meta data in redis database on compute nodes')
    parser.add_argument('configfile', type=str,
            help='Config file used to program snaps')
    parser.add_argument('-s', dest='snaphosts', nargs='+', type=str,
            required=True, help='fpga host names')
    parser.add_argument('-r', dest='dry_run', action='store_true',
            help='Dry run (do not publish)')
    parser.add_argument('-i', dest='redis_instance',
            help='specify the instance enumeration of the redis gateway',
            default=0, type=int)
    parser.add_argument('-I', dest='ignore_control',
            help='ignore the sky frequency from the antennas',
            action='store_true')

    args = parser.parse_args()

    r = redis.Redis(host=REDISHOST)

    with open(args.configfile, 'r') as fh:
        config = yaml.load(fh, Loader=yaml.SafeLoader)

    dest_port = config.get('dest_port')

    voltage_config = config.get('voltage_output')

    dests      = voltage_config['dests']
    n_chans    = voltage_config['n_chans']
    start_chan = voltage_config['start_chan']
    nants      = len(args.snaphosts)
    n_dests    = len(dests)
    sync_time  = int(r.get('SYNCTIME'))
    snapseq    = ",".join([isnap.replace('frb-snap','').replace('-pi','')
        for isnap in args.snaphosts]) #this contains the "physical" snapID

    n_chans_per_dest = n_chans // n_dests

    mapping = get_channel_selection(dests, start_chan,
            n_chans_per_dest)

    skyfreq_mapping, antname_mapping = get_snap_mapping(args.snaphosts,
            args.ignore_control)

    if len(skyfreq_mapping.values()) != 1:
        sys.stderr.write("WARNING: antennas are tuned to different freqs, "
                "OBSFREQ will be given the first value of OBSNFREQ\n")

    obsfreq = skyfreq_mapping[args.snaphosts[0]]

    for ip,lst in mapping.items():
        n_packets_per_dest = int(np.ceil(n_chans_per_dest / MAX_CHANS_PER_PKT))
        n_chans_per_pkt  = n_chans_per_dest // n_packets_per_dest
        schan = lst[0]
        nstrm = n_chans_per_dest // n_chans_per_pkt

        chan_bw = FOFF
        obsbw   = len(lst)*chan_bw

        host = socket.gethostbyaddr(ip)[0].replace("-40g","")
        channel_name = REDISSETGW.substitute(host=host, inst=args.redis_instance)

        # these are instance specific
        key_val = {
                'OBSBW'    : obsbw,
                'SCHAN'    : schan,
                'NSTRM'    : nstrm,
                'OBSFREQ'  : obsfreq,
                'BINDHOST' : BINDHOST,		# cannot change once hashpipe has started.
                'BINDPORT' : dest_port,		# cannot change once hashpipe has started.
								'BANK'     : BANK,
								'SNAPSEQ'  : snapseq,
        }

        key_val_str = "\n".join(['%s=%s' %(key,val)
            for key,val in key_val.items()])
        print("channel_name:", channel_name)
        print(key_val_str)
        print()

        if not args.dry_run:
            # publish them values
            r.publish(channel_name, key_val_str)

    # these are common for all nodes
    key_val = {
            'CHAN_BW'  : chan_bw,
            'FENCHAN'  : FENCHAN,
            'NANTS'    : nants,
            'NPOL'     : NPOL,
            'PKTNCHAN' : n_chans_per_pkt,
            'TBIN'     : TBIN,
            'NBITS'    : NBITS,
            'PKTNTIME' : N_TIMES_PER_PKT,
            'SYNCTIME' : sync_time,
            'PROJID'   : PROJID,
            'BACKEND'  : BACKEND,
            'DATADIR'  : DATADIR,
            'PKTFMT'   : PKTFMT,
            'SNAPPAT'  : SNAPPAT,
            'ANTNAMES' : '1a,1f'
    }

    key_val_str = "\n".join(['%s=%s' %(key,val)
        for key,val in key_val.items()])

    print("channel_name:", channel_name)
    print(key_val_str)
    if not args.dry_run:
        # publish them values
        r.publish(REDISSET, key_val_str)


if __name__ == "__main__":
    main()
