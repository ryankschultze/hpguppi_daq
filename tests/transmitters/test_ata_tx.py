#!/usr/bin/env python

import socket
import time
import struct
import argparse
import numpy as np

NPOL = 2
NCHAN_MAX = 256
NTIME = 16
NBIT = 8 #real+imag

parser = argparse.ArgumentParser(description='Emulate a bunch of ATA SNAPs sending data',
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-p', '--port', type=int, default=10000,
                    help='UDP port to which data should be sent')
parser.add_argument('-i', '--ip', type=str, default='100.100.100.100',
                    help='IP address to which data should be sent')
parser.add_argument('-f', '--nfeng', type=int, default=10,
                    help='Number of F-engines to emulate')
parser.add_argument('-c', '--nchan', type=int, default=256,
                    help='Number of frequency channels to send')
parser.add_argument('-g', '--gbps', type=float, default=None,
                    help='Target number of Gbps. Will throttle if necessary. If None, send at max speed')
args = parser.parse_args()

sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

# Data have format (slowest to fastest): channel x time x polarization x complexity x 4 bits

nchan_block = (args.nchan + NCHAN_MAX - 1) // NCHAN_MAX
if args.nchan > NCHAN_MAX:
    if args.nchan % nchan_block:
        print("ERROR: Trying to send %d channels in %d blocks, which doesn't divide evenly." % (args.nchan, nchan_block))
        print("Maximum number of channels per block is %d" % NCHAN_MAX)
        exit()

nchan_per_block = args.nchan // nchan_block

print("Sending %d channels from %d F-Engines" % (args.nchan, args.nfeng))
print("Sending as %d blocks of %d channels" % (nchan_block, nchan_per_block))

payload_len = nchan_per_block * NTIME * NPOL * NBIT // 8
print("Payload length is %d bytes" % payload_len)
# use zeros as the payload for now
payload = b"\x00" * payload_len

# Packet format
#   struct packet{
#     uint8_t version;
#     uint8_t type;
#     uint16_t nchan;
#     uint16_t chan;
#     uint16_t feng_id;
#     uint64_t timestamp;
#     complex4 data[nchan_per_block, NTIME, NPOL]; // 4-bit real + 4-bit imaginary
#   };

# Header fields which never change
h_version = 0
h_type = 0
h_nchan = nchan_per_block

pkt_count = 0
t = 0
tick = time.time()
throttle_tick = time.time()
NPKT_THROTTLE = 100
nbytes_throttle = NPKT_THROTTLE * payload_len
npkt_sent = 0
while(True):
    for f in range(args.nfeng):
        for c in range(nchan_block):
            header = struct.pack(">BBHHHQ", h_version, h_type, h_nchan, c*nchan_per_block, f, t)
            sock.sendto(header + payload, (args.ip, args.port))
            pkt_count += 1
            if args.gbps is None:
                continue
            if pkt_count % NPKT_THROTTLE == 0:
                    throttle_tock = time.time()
                    elapsed = throttle_tock - throttle_tick
                    target_elapsed = 8*nbytes_throttle / 1e9 / args.gbps
                    if target_elapsed > elapsed:
                        time.sleep(target_elapsed - elapsed)
                    throttle_tick = time.time()

    t += NTIME

    if time.time() - tick > 1:
        tock = time.time()
        elapsed = tock-tick
        npkt_sent_new = pkt_count - npkt_sent
        nbyte_sent_new = npkt_sent_new * payload_len
        print("Sent %d packets (%.1f MB) in %.2fs (%.2fGb/s)" % (npkt_sent_new, nbyte_sent_new / 1.0e6, elapsed, 8*nbyte_sent_new / elapsed / 1e9))
        tick = time.time()
        npkt_sent = pkt_count


