#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import sys
from collections import defaultdict
from Bio import SeqIO
import pysam
import argparse
from pathlib import Path
import csv
csv.field_size_limit(1000000000)


def main():
    parser = argparse.ArgumentParser(description='KIR reference trimmer')
    parser.add_argument('sam', type=Path, help='KIR SAM filename')
    parser.add_argument('--name', type=str,
                        help='primary sequence name', required=True)
    parser.add_argument('--start', type=int, required=True,
                        help='start position (0-origin) in primary')
    parser.add_argument('--end', type=int, required=True,
                        help='end position (0-origin) in primary')
    args = parser.parse_args()

    # breakpoints[qname] = ([start], [end])
    breakpoints = defaultdict(lambda: (set(), set()))

    sam = pysam.AlignmentFile(args.sam, "r")
    for x in sam:
        tstart = x.reference_start
        tend = x.reference_end
        qstart = x.query_alignment_start
        qend = x.query_alignment_end
        is_reverse = x.is_reverse
        qname = x.query_name
        tname = x.reference_name
        if tname == args.name:
            if tstart <= args.start <= tend:
                qbreakpoint = next(
                    (qpos for (qpos, tpos) in x.get_aligned_pairs() if tpos == args.start), None)
                print(qname, 'start', tstart, tend, qstart,
                      qend, is_reverse, qbreakpoint, file=sys.stderr)
                breakpoints[qname][0].add(qbreakpoint)
            if tstart <= args.end <= tend:
                # print('end', qname)
                qbreakpoint = next(
                    (qpos for (qpos, tpos) in x.get_aligned_pairs() if tpos == args.end), None)
                print(qname, 'end', tstart, tend, qstart,
                      qend, is_reverse, qbreakpoint, file=sys.stderr)
                breakpoints[qname][1].add(qbreakpoint)

    for (qname, (starts, ends)) in breakpoints.items():
        print(qname, list(starts), list(ends), sep='\t')


if __name__ == '__main__':
    main()
