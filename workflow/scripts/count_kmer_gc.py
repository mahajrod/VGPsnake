#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import bz2
import sys
import gzip
import argparse

if sys.version_info[0] == 3:
    from io import TextIOWrapper as file


def metaopen(filename, flags, buffering=None, compresslevel=5):
    if not isinstance(filename, str): # or isinstance(filename, gzip.GzipFile) or isinstance(filename, bz2.BZ2File):
        if isinstance(filename, file):
            return filename
        else:
            raise ValueError("ERROR!!! Not file object or str: {}".format(str(filename)))
    elif filename[-3:] == ".gz":
        return gzip.open(filename, flags + ("t" if "b" not in flags else ""), compresslevel=compresslevel)
    elif filename[-4:] == ".bz2":
        return bz2.open(filename, flags + ("t" if "b" not in flags else ""), compresslevel=compresslevel)
    else:
        if buffering is not None:
            return open(filename, flags, buffering=buffering)
        else:
            return open(filename, flags)

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input 2-column (kmer and count) tab-separated file with kmer counts. Default: stdin")
parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output 2-column (gc content and count) tab-separated file with. Default: stdout")

args = parser.parse_args()

with metaopen(args.input, "r", buffering=100000000) as in_fd, metaopen(args.output, "w", buffering=100000000) as out_fd:
    for line in in_fd:
        line_list = line.split() # splits on any space-character and trims \n from thew end
        out_fd.write("{0}\t{1}\n".format(line_list[0].count("G") + line_list[0].count("C"), line_list[1]))

