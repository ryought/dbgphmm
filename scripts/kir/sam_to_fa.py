#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
from Bio import SeqIO
import pysam
import argparse
from pathlib import Path
import csv
csv.field_size_limit(1000000000)


def main():
    parser = argparse.ArgumentParser(description='Read trimmer')
    parser.add_argument('sam', type=Path, help='Read-vs-ref SAM filename')
    parser.add_argument('fasta', type=Path, help='Read FASTA')
    parser.add_argument('--target', type=str,
                        help='target sequence (reference chromosome) name')
    parser.add_argument('--start', type=int, help='start position (0-origin)')
    parser.add_argument('--end', type=int, help='end position (0-origin)')
    args = parser.parse_args()

    sam = pysam.AlignmentFile(args.sam, "r")

    fasta = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

    for x in sam.fetch(args.target, args.start, args.end):
        tstart = x.reference_start
        tend = x.reference_end
        qstart = x.query_alignment_start
        qend = x.query_alignment_end
        is_reverse = x.is_reverse
        qname = x.query_name
        query = fasta[qname].seq
        qlen = len(query)

        if x.is_secondary or x.is_supplementary or x.is_unmapped:
            continue

        if x.is_reverse:
            continue

        if tstart <= args.start <= tend:
            # case1: have start breakpoint
            qbreakpoint = next(
                (i for (i, j) in x.get_aligned_pairs() if j == args.start), None)
            # print(x)
            # print(tstart - args.start, tend - args.start, qstart, qend, is_reverse, qbreakpoint)
            print('>{}:{}-{}:+:start'.format(qname, qbreakpoint, qlen))
            print(query[qbreakpoint:])
        elif tstart <= args.end <= tend:
            # case2: have end breakpoint
            qbreakpoint = next(
                (i for (i, j) in x.get_aligned_pairs() if j == args.end), None)
            print('>{}:{}-{}:+:end'.format(qname, 0, qbreakpoint))
            print(query[:qbreakpoint])
        else:
            # case3: no breakpoint
            print('>{}:{}-{}:+:all'.format(qname, 0, qlen))
            print(query)


if __name__ == '__main__':
    main()
