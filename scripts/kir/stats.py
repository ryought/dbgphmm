#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
from pafpy import PafFile
import argparse
from pathlib import Path
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('paf', type=Path, help='PAF')
    args = parser.parse_args()

    # precision = n-nm / n
    # where
    # n = total assembled base
    # nm = total mismatch base
    nms = dict()

    # recall = a / g
    # g = total bases of genome
    # a = total bases of genome which has alignment of reads
    # coverage[hap] = [0,0,0,0....]
    coverage = dict()

    with PafFile(args.paf) as paf:
        for record in paf:
            nm = record.get_tag("NM").value
            # print(record)
            # print(record.qname, record.qlen, nm)

            if record.tname not in coverage:
                coverage[record.tname] = [0 for _ in range(record.tlen)]
            for i in range(record.tstart, record.tend):
                coverage[record.tname][i] = 1

            nms[record.qname] = (nm, record.qlen)

    for hap in coverage.keys():
        n = sum(coverage[hap])
        m = len(coverage[hap])
        print(hap, n, m, n / m)
    for contig in nms.keys():
        n, m = nms[contig]
        print(contig, n, m, (m - n) / m)


if __name__ == '__main__':
    main()
