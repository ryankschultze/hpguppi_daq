#!/usr/bin/env python
import subprocess
import argparse
import time
import glob
import os
import re

from hashpipe_aux import *

parser = argparse.ArgumentParser(description='Starts up a hashpipe, emulated localhost transmission '
                                             'and verification of the RAW capture.',
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# hashpipe status arguments
parser.add_argument('-P', '--projid', type=str, default='tx_test',
                    help='The PROJID hashpipe status value')
parser.add_argument('-I', '--instance', type=int, default=0,
                    help='Instance ID of the hashpipe.')
# transmission arguments
parser.add_argument('-p', '--port', type=int, default=10000,
                    help='UDP port to which data should be sent')
parser.add_argument('-i', '--ip', type=str, default='100.100.100.100',
                    help='IP address to which data should be sent')
parser.add_argument('-f', '--nfeng', type=int, default=10,
                    help='Number of F-engines to emulate')
parser.add_argument('-c', '--nchan', type=int, default=256,
                    help='Number of frequency channels to send')
parser.add_argument('-t', '--starttime', type=int, default=0,
                    help='Timestamp of first packet sent')
parser.add_argument('-C', '--startchan', type=int, default=0,
                    help='Start channel of first packet sent')
parser.add_argument('--ntime', type=int, default=None, required=True,
                    help='If set, generate (or send) `ntime` time samples of data')
parser.add_argument('--prepackets', type=int, default=1000,
                    help='Number of packets to send before starttime')
parser.add_argument('--postpackets', type=int, default=1000,
                    help='Number of packets to send after starttime+ntime')
parser.add_argument('-g', '--gbps', type=float, default=None,
                    help='Target number of Gbps. Will throttle if necessary. If None, send at max speed')
parser.add_argument('--infile', type=str, default=None, required=True,
                    help='File containing data to be sent. Data should be in time x feng x chan x complexity order')
parser.add_argument('-r', '--reusehashpipe', action='store_true',
                    help='Target number of Gbps. Will throttle if necessary. If None, send at max speed')
args = parser.parse_args()

NCHAN_MAX = 256
nstrm = (args.nchan + NCHAN_MAX - 1) // NCHAN_MAX
if args.nchan > NCHAN_MAX:
    if args.nchan % nstrm:
        print('ERROR: Trying to send %d channels in %d streams, which doesn\'t divide evenly.' % (args.nchan, nstrm))
        print('Maximum number of channels per block is %d' % NCHAN_MAX)
        exit()
nchan_per_pkt = args.nchan // nstrm
# --ntime 4000000 -c 1440 -f 3 --infile noise_4Mt_1440c_3f.dat -i 10.11.1.155 -t 1978750000 --startchan 2048

# print("\nKilling existing hashpipes")
# subprocess.run(['sudo', 'pkill', '-e', '-f', '.*/hashpipe\s.*'])


# subprocess.run(['hashpipe_check_status', '--key=PKTNCHAN', '--int='+str(256)]) 
# subprocess.run(['hashpipe_check_status', '--key=NANTS',    '--int='+str(1)]) 
# subprocess.run(['hashpipe_check_status', '--key=NSTRM',    '--int='+str(1)]) 
# subprocess.run(['hashpipe_check_status', '--key=OBSSTART', '--int='+str(0)]) 
# subprocess.run(['hashpipe_check_status', '--key=OBSSTOP',  '--int='+str(0)]) 

for hashpipe_start_count in range(2):
    if not args.reusehashpipe:
        ##### Start Hashpipe ##### 
        print('\n######Starting Hashpipe#####\n')
        start_hashpipe(args.instance, 'lo') # assume it will kill existing hashpipes
        time.sleep(2)

    ##### Start Hashpipe Gateway #####
    start_redis_gateway(args.instance)
    time.sleep(3)

    if not block_until_pulse_change(20):
        kill_hashpipe_relevant()
        clear_shared_memory()
        if hashpipe_start_count > 1:
            exit()
    else:
        time.sleep(5)
        break

##### Set Hashpipe Keys #####
print('\n######Setting Hashpipe Status Keys#####\n')

subprocess.run(['hashpipe_check_status', '--key=PKTFMT',  '--string=ATASNAPV']) 
subprocess.run(['hashpipe_check_status', '--key=BACKEND',  '--string=tx_test']) 
subprocess.run(['hashpipe_check_status', '--key=BANKNAM',  '--int=0']) 
subprocess.run(['hashpipe_check_status', '--key=NBITS',    '--int=4']) 
subprocess.run(['hashpipe_check_status', '--key=NPOL',     '--int=2']) 
subprocess.run(['hashpipe_check_status', '--key=PKTNTIME', '--int=16'])
subprocess.run(['hashpipe_check_status', '--key=CHAN_BW',  '--float=0.25']) 
subprocess.run(['hashpipe_check_status', '--key=OBSBW',    '--double=128.0']) 
subprocess.run(['hashpipe_check_status', '--key=OBSFREQ',  '--int=1420']) 
# subprocess.run(['hashpipe_check_status', '--key=OBSNCHAN',  '--int='+str(nchan_per_pkt*args.nfeng*nstrm)]) 

subprocess.run(['hashpipe_check_status', '--key=PROJID',   '--string='+args.projid])
subprocess.run(['hashpipe_check_status', '--key=PKTNCHAN', '--int='+str(nchan_per_pkt)]) 
subprocess.run(['hashpipe_check_status', '--key=NANTS',    '--int='+str(args.nfeng)]) 
subprocess.run(['hashpipe_check_status', '--key=NSTRM',    '--int='+str(nstrm)]) 
subprocess.run(['hashpipe_check_status', '--key=OBSSTART', '--int='+str(args.starttime)]) 
subprocess.run(['hashpipe_check_status', '--key=OBSSTOP',  '--int='+str(args.starttime + args.ntime)]) 
subprocess.run(['hashpipe_check_status', '--key=SCHAN',    '--int='+str(args.startchan)]) 
# subprocess.run(['hashpipe_check_status', '--key=TBIN',     '--double='+str(nchan_per_pkt*args.nfeng*nstrm / 128.0 / 1e6)]) # 6.3164065f-5


##### Transmit Test Vector #####
transmissionCmd = ['numactl', '--cpunodebind=0', '--membind=0', 
                   './transmitters/test_ata_tx.py'
                ]
transmissionCmd.append('--port')
transmissionCmd.append(str(args.port))
transmissionCmd.append('--ip')
transmissionCmd.append(args.ip)
transmissionCmd.append('--nfeng')
transmissionCmd.append(str(args.nfeng))
transmissionCmd.append('--nchan')
transmissionCmd.append(str(args.nchan))
transmissionCmd.append('--starttime')
transmissionCmd.append(str(args.starttime))
transmissionCmd.append('--prepackets')
transmissionCmd.append(str(args.prepackets))
transmissionCmd.append('--postpackets')
transmissionCmd.append(str(args.postpackets))
transmissionCmd.append('--startchan')
transmissionCmd.append(str(args.startchan))
transmissionCmd.append('--infile')
transmissionCmd.append(args.infile)
transmissionCmd.append('--ntime')
transmissionCmd.append(str(args.ntime))
if args.gbps:
    transmissionCmd.append('--gbps')
    transmissionCmd.append(str(args.gbps))

print('\n######Transmitting test vector#####\n')
time.sleep(5)
print(transmissionCmd)
subprocess.run(transmissionCmd) 

##### Verify Hashpipe Capture #####

firstrawFilepath = get_latest_raw_stem_in_dir(get_hashpipe_capture_dir())+'.0000.raw'
rawtime = time.ctime(os.path.getmtime(firstrawFilepath))
print(firstrawFilepath, '\nLast modification time:', rawtime, '(Local time)\n')


verificationCmd = ['./verification/verify_ata.py']
verificationCmd.append('--rawfile')
verificationCmd.append(firstrawFilepath)
verificationCmd.append('--infile')
verificationCmd.append(args.infile)
verificationCmd.append('--nfeng')
verificationCmd.append(str(args.nfeng))
verificationCmd.append('--nchan')
verificationCmd.append(str(args.nchan))

print(' '.join(verificationCmd))
subprocess.run(verificationCmd)

print(firstrawFilepath)