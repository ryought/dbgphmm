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
    parser = argparse.ArgumentParser(description='KIR chunk lifter')
    parser.add_argument(
        'sam', type=Path, help='KIR SAM filename (map chunk onto sequences)')
    parser.add_argument('--pos', type=int, required=True,
                        help='position (0-origin) in chunk')
    args = parser.parse_args()

    # breakpoints[qname] = ([start], [end])
    # breakpoints = defaultdict(lambda: (set(), set()))

    sam = pysam.AlignmentFile(args.sam, "r")
    for x in sam:
        tstart = x.reference_start
        tend = x.reference_end
        qstart = x.query_alignment_start
        qend = x.query_alignment_end
        is_reverse = x.is_reverse
        qname = x.query_name
        tname = x.reference_name
        if qstart <= args.pos <= qend:
            tbreakpoint = next(
                (tpos for (qpos, tpos) in x.get_aligned_pairs() if qpos == args.pos), None)
            print(tname, tbreakpoint, sep='\t')


if __name__ == '__main__':
    main()
