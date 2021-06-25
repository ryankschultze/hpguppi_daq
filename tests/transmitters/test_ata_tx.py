#!/usr/bin/env python

import os
import socket
import time
import struct
import argparse
import numpy as np

# Packet header flags
PKT_IS_VOLTAGE = 1<<7 # packet contains voltage data
PKT_VER = 2 << 5
PKT_TYPE = 0


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
parser.add_argument('-t', '--starttime', type=int, default=0,
                    help='Timestamp of first packet sent')
parser.add_argument('-C', '--startchan', type=int, default=0,
                    help='Start channel of first packet sent')
parser.add_argument('--ntime', type=int, default=None,
                    help='If set, generate (or send) `ntime` time samples of data')
parser.add_argument('--outfile', type=str, default=None,
                    help='If set, generate a data file with T samples and then exit')
parser.add_argument('-g', '--gbps', type=float, default=None,
                    help='Target number of Gbps. Will throttle if necessary. If None, send at max speed')
parser.add_argument('--infile', type=str, default=None,
                    help='File containing data to be sent. Data should be in time x feng x chan x complexity order')
parser.add_argument('-m', '--miss', type=int, default=None,
                    help='If set, specify a number of packets to be sent, after which one will be deliberately dropped')
args = parser.parse_args()

assert NBIT == 8, "Only 4+4 bit data currently supported"

if args.outfile is not None:
    print("Not sending data. Generating a test file instead")
    if args.ntime is None:
        print("ERROR: When generating a test file the -T option must be specified")
        exit()
    if not (args.ntime % NTIME == 0):
        print("ERROR: The number of time samples generated (set with -T) must be divisible by %d" % NTIME)
        exit()
    shape = [args.ntime // NTIME, args.nfeng, args.nchan, NTIME, NPOL]
    nbytes = np.prod(shape)
    print("Generating %d bytes random data with axes [Time Blocks x F-engines x Channels x Time x Pols]: %s" % (nbytes, shape))
    d = np.random.randint(0, 255, size=shape, dtype=np.uint8)
    with open(args.outfile, "wb") as fh:
        fh.write(d.tobytes())
    exit()

sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

# Data have format (slowest to fastest): channel x time x polarization x complexity x 4 bits

nchan_block = (args.nchan + NCHAN_MAX - 1) // NCHAN_MAX
if args.nchan > NCHAN_MAX:
    if args.nchan % nchan_block:
        print("ERROR: Trying to send %d channels in %d blocks, which doesn't divide evenly." % (args.nchan, nchan_block))
        print("Maximum number of channels per block is %d" % NCHAN_MAX)
        exit()

nchan_per_block = args.nchan // nchan_block

if args.infile is not None:
    if not os.path.isfile(args.infile):
        print("ERROR: file %s does not exist" % args.infile)
        exit()
    print("Reading file %s" % args.infile)
    with open(args.infile, "rb") as fh:
        d_raw = fh.read()
    nbyte = len(d_raw)
    assert nbyte % (args.nfeng * args.nchan * NTIME * NPOL * NBIT // 8) == 0, "File does not contain an integer number of packets"
    ntime = nbyte // (args.nfeng * args.nchan * NPOL * NBIT // 8)
    print("Read %d bytes (%d times)" % (nbyte, ntime))
    ntime_block = ntime // NTIME
    d = np.frombuffer(d_raw, dtype=np.uint8).reshape([ntime_block, args.nfeng, nchan_block, nchan_per_block, NTIME, NPOL])


print("Sending %d channels from %d F-Engines" % (args.nchan, args.nfeng))
print("Sending as %d blocks of %d channels" % (nchan_block, nchan_per_block))

payload_len = nchan_per_block * NTIME * NPOL * NBIT // 8
print("Payload length is %d bytes" % payload_len)
# use zeros as the payload for now
payload = b"\x00" * payload_len

# Packet format
#   struct packet{
#     uint8_t version;    // Packet version. MSB is always 1 for voltage packets
#     uint8_t type;       // Type. Currently 0 (indicating channel is slowest axis, 4+4 bit data, no GPU swizzle)
#     uint16_t nchan;     // Number of channels in this packet
#     uint16_t chan;      // First channel in this packet
#     uint16_t feng_id;   // 0-indexed F-engine ID
#     uint64_t timestamp; // Timestamp of first spectra in this packet
#     complex4 data[nchan_per_block, NTIME, NPOL]; // 4-bit real + 4-bit imaginary
#   };

# Header fields which never change
h_version = PKT_IS_VOLTAGE + PKT_VER
h_type = PKT_TYPE
h_nchan = nchan_per_block

pkt_count = 0
t = 0
tn = 0 # Count of complete sets of all fengs / chans
tick = time.time()
throttle_tick = time.time()
NPKT_THROTTLE = 100
nbytes_throttle = NPKT_THROTTLE * payload_len
npkt_sent = 0
sending_0s = True
print('Sending payload of 0s')
while(True):
    for f in range(args.nfeng):
        for c in range(nchan_block):
            if args.miss is not None:
                if pkt_count % args.miss == 0:
                    #skip this packet
                    pkt_count += 1
                    continue
            header = struct.pack(">BBHHHQ", h_version, h_type, h_nchan, args.startchan + c*nchan_per_block, f, t)
            if args.infile and t >= args.starttime and tn < ntime_block:
                if sending_0s:
                    print('Sending payload of {} (t={})'.format(args.infile, t))
                sending_0s = False
                payload = d[tn, f, c, :, :, :].tobytes()
            # else: send payload of 0s
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
    if not sending_0s:
        tn += 1
    if args.ntime is not None:
        if t > args.ntime:
            print("Exiting because %d time samples have been sent" % (t))
            exit()

    if args.infile is not None:
        if tn >= ntime_block and not sending_0s:
            print("Sending zero payloads because all %d time samples in the input file have been sent (t=%d)" % (tn*NTIME, t))
            payload = b"\x00" * payload_len
            sending_0s = True

    if time.time() - tick > 1:
        tock = time.time()
        elapsed = tock-tick
        npkt_sent_new = pkt_count - npkt_sent
        nbyte_sent_new = npkt_sent_new * payload_len
        print("Sent %d packets (%.1f MB) in %.2fs (%.2fGb/s)" % (npkt_sent_new, nbyte_sent_new / 1.0e6, elapsed, 8*nbyte_sent_new / elapsed / 1e9))
        tick = time.time()
        npkt_sent = pkt_count


