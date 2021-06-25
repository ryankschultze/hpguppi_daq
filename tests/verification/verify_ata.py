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


HEADER_KEY_VAL_SIZE = 80 #bytes
DIRECT_IO_SIZE      = 512

class Guppi():
    """
    A very basic guppi raw file reader
    """
    def __init__(self, fname):
        if type(fname) != str:
            raise RuntimeError("Please provide string filename")
        self.fname = fname
        self.file  = open(fname, "rb")

    def __del__(self):
        self.file.close()

    #@jit(nopython=True)
    def _parse_header(self):
        header = {}
        nbytes_read = 0
        hread = self.file.read(HEADER_KEY_VAL_SIZE).decode('UTF-8')
        if not hread: # we have reachec the end of file
            return None
        nbytes_read += HEADER_KEY_VAL_SIZE
        while not hread.startswith("END"):
            key, val = hread.split("=")
            key = key.strip()
            val = val.strip()

            try:
                if "." in val:
                    val = float(val)
                else:
                    val = int(val)
            except ValueError:
                val = val.strip("'").strip()

            header[key] = val
            hread = self.file.read(HEADER_KEY_VAL_SIZE).decode('UTF-8')
            nbytes_read += HEADER_KEY_VAL_SIZE

        assert hread == "END"+" "*77, "Not a GUPPI RAW format"
        if header['DIRECTIO']:
            remainder = nbytes_read % DIRECT_IO_SIZE
            to_seek = (DIRECT_IO_SIZE - remainder)%DIRECT_IO_SIZE
            _ = self.file.read(to_seek)
            nbytes_read += to_seek
        header['HEADER_SIZE'] = nbytes_read
        return header

    #@jit(nopython=True)
    def read_next_block(self, return_8bit=False):
        header = self._parse_header()
        if not header:
            return None, None

        npol     = header['NPOL']
        obsnchan = header['OBSNCHAN']
        nbits    = header['NBITS']
        blocsize = header['BLOCSIZE']
        try:
            nants    = header['NANTS']
        except KeyError as e:
            nants = -1

        if nbits != 4:
            raise NotImplementedError("Only 4-bit data is implemented")

        if return_8bit:
            data = np.fromfile(self.file, dtype=np.uint8, count=blocsize)
        else:
            data_raw = np.fromfile(self.file, dtype=np.int8, count=blocsize)
            data = np.zeros_like(data_raw, dtype=np.complex64)
            data[:] = (data_raw >> 4) + 1j*(data_raw << 4 >> 4)

        nsamps_per_block = int(blocsize / (2*npol * obsnchan * (nbits/8.0)))

        if (2 * npol * obsnchan * (nbits/8.0) * nsamps_per_block ) != blocsize:
            raise RuntimeError("Bad block geometry: 2*%i*%i*%f*%i != %i"\
                    %(npol, obsnchan, nbits/8.0, nsamps_per_block, blocsize))

        if nants != -1: # "multi-antenna" raw file
            nchan_per_ant    = obsnchan//nants

            if (nchan_per_ant * nants) != obsnchan:
                raise RuntimeError("obsnchan does not equally divide across antennas: "\
                        "obsnchan: %i, nants: %i, nchan_per_ant: %i",
                        obsnchan, nants, nchan_per_ant)

            data = data.reshape(nants, nchan_per_ant, nsamps_per_block, npol)
        else:
            data = data.reshape(obsnchan, nsamps_per_block, npol)

        if header['DIRECTIO']:
            remainder = blocsize % DIRECT_IO_SIZE
            to_seek = (DIRECT_IO_SIZE - remainder)%DIRECT_IO_SIZE
            if to_seek:
                _ = self.file.read(to_seek)

        return header, data

parser = argparse.ArgumentParser(description='Emulate a bunch of ATA SNAPs sending data',
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-f', '--nfeng', type=int, default=10,
                    help='Number of F-engines to emulate')
parser.add_argument('-c', '--nchan', type=int, default=256,
                    help='Number of frequency channels to send')
parser.add_argument('--rawfile', type=str, default=None, required=True,
                    help='.raw file recorded by pipeline')
parser.add_argument('--ntime', type=int, default=None,
                    help='Upper limit of `ntime` time samples of data to verify')
parser.add_argument('--infile', type=str, default=None, required=True,
                    help='File containing data which was sent. Data should be in time x feng x chan x complexity order')
args = parser.parse_args()

assert NBIT == 8, "Only 4+4 bit data currently supported"


if not os.path.isfile(args.infile):
    print("ERROR: file %s does not exist" % args.infile)
    exit()
if not os.path.isfile(args.rawfile):
    print("ERROR: file %s does not exist" % args.rawfile)
    exit()

# doesn't seem like we can read sequentially as we verify, as the raw data isn't in the same order
print("Reading file %s" % args.infile)
with open(args.infile, "rb") as fh:
    #d_raw = fh.read(10000 * (args.nfeng * args.nchan * NTIME * NPOL * NBIT // 8))
    d_raw = fh.read()
nbyte = len(d_raw)
assert nbyte % (args.nfeng * args.nchan * NTIME * NPOL * NBIT // 8) == 0, "File does not contain an integer number of packets"
ntime = nbyte // (args.nfeng * args.nchan * NPOL * NBIT // 8)
print("Read %d bytes (%d times)" % (nbyte, ntime))
assert ntime % NTIME == 0
ntime_block = ntime // NTIME
d = np.frombuffer(d_raw, dtype=np.uint8).reshape([ntime_block, args.nfeng, args.nchan, NTIME, NPOL])

bad_cnt = 0
good_cnt = 0
g = Guppi(args.rawfile)
n_block = 0
n_bytes_checked = 0
capture_offset = -1
while(True):
    print("Testing block %d (checked %.1f Mbytes)" % (n_block, n_bytes_checked / 1.0e6))
    gh, gd = g.read_next_block(return_8bit=True)
    if gh is None:
        break
    #if n_block > 100:
    #    break
    g_nants, g_nchans, g_ntimes, g_npols = gd.shape
    print('gd.shape', gd.shape)
    print("PKTIDX:", gh['PKTIDX'])
    if capture_offset == -1:
        capture_offset = gh['PKTIDX'] - gh['PKTSTART']
        print('capture_offset:', capture_offset)
    # gd has shape: nants x nchans x ntimes x npols
    # infile data has shape: ntime x nants x nchans x NTIME x NPOL
    assert g_ntimes % NTIME == 0
    # TODO: packets up until the first whose index is a block boundary (N*PIPERBLK)
    # are not captured. ensure that the transmission is block aligned, or capture the 
    # first PKTIDX of the transmission somehow so that we can ignore the relevant number
    # of samples when verifying
    block_time = (n_block) * g_ntimes
    for t in range(g_ntimes):
        golden_slice = gd[:, :, t, :]
        record_slice = d[(block_time + t + capture_offset)//NTIME, :, :, (block_time + t + capture_offset)%NTIME, :]
        diff = np.equal(golden_slice, record_slice)
        if diff.all():
            good_cnt += g_nants * g_nchans * g_npols
        else:
            # print('-'*20+'\n', n_block, t, '\n', golden_slice, '\nvs\n', record_slice, '\n'+'-'*20)
            print('-'*20+'\n', n_block, t, diff, '\n'+'-'*20)
            bad_cnt += g_nants * g_nchans * g_npols
            break
        # TODO why are samples after 12239 bad?
        #if t == 12239:
        #    break
    n_block += 1
    n_bytes_checked += g_nants * g_nchans * g_ntimes * g_npols
    print("Good count: %d; Bad count %d;" % (good_cnt, bad_cnt))
    if bad_cnt != 0:
        break

print('capture_offset:', capture_offset)
if bad_cnt == 0:
    print("TEST PASSED")
else:
    print("TEST FAILED")

