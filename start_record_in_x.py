#!/home/obsuser/miniconda3/envs/ATAobs/bin/python
import numpy as np
import argparse
import redis
import time
import socket
from string import Template

REDISHOST = 'redishost'
REDISSET  = 'hashpipe:///set'

FENCHAN   = 4096
BW        = 1.024e9
TBIN      = 1/(BW/FENCHAN)

DEFAULT_START_IN = 2
DEFAULT_OBS_TIME = 300.


def main():
    parser = argparse.ArgumentParser(description='start observation '
            'in x seconds')
    parser.add_argument('-n', type=float, default=DEFAULT_OBS_TIME,
            help='Total obs time [%.1f]' %DEFAULT_OBS_TIME)
    parser.add_argument('-i', type=int, default=DEFAULT_START_IN,
            help='Seconds from now to start obs [%i]' %DEFAULT_START_IN)
    parser.add_argument('-r', action='store_true',
            help='Reset OBSSTART and OBSSTOP to 0')
    parser.add_argument('-d', action='store_true',
            help='dry run (don\'t publish)')
    parser.add_argument('-H', type=int, default=-1,
            help='hashpipe instance id to target [-1: for broadcast]')

    args = parser.parse_args()

    r = redis.Redis(host=REDISHOST)

    sync_time = int(r.get("SYNCTIME"))

    if args.r:
        cmd = "OBSSTART=0\nOBSSTOP=0"
    else:
        t_now  = time.time()
        t_in_2 = int(np.ceil(t_now + args.i))

        tdiff = t_in_2 - sync_time

        obsstart = int(tdiff/TBIN)

        npckts_to_record = int(args.n/TBIN)
        obsstop = obsstart + npckts_to_record
        print("OBSSTART: %i" %obsstart)
        print("OBSSTOP : %i" %obsstop)
        cmd = "OBSSTART=%i\nOBSSTOP=%i"  %(obsstart, obsstop)

    if args.H > -1:
        host = socket.gethostname()
        REDISSETGW = Template('hashpipe://${host}/${inst}/set')
        REDISSET = REDISSETGW.substitute(host=host, inst=args.H)
        print(REDISSET)


    if args.d:
        print (cmd)
    else:
        r.publish(REDISSET, cmd)


if __name__ == "__main__":
    main()
