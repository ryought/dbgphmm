#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import csv
csv.field_size_limit(1000000000)


def main():
    parser = argparse.ArgumentParser(description='Dotplot from paf')
    parser.add_argument('paf', type=Path, help='PAF filename')
    parser.add_argument('--query', type=str, help='query sequence name')
    parser.add_argument('--target', type=str,
                        help='target (reference) sequence name')
    parser.add_argument('--qs', type=int, help='query start')
    parser.add_argument('--qe', type=int, help='query end')
    parser.add_argument('--ts', type=int, help='target start')
    parser.add_argument('--te', type=int, help='target end')

    parser.add_argument('--out', type=str, help='plot save destination')
    args = parser.parse_args()

    plt.axes().set_aspect('equal')
    with open(args.paf) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0] == args.query and row[5] == args.target:
                qlen = int(row[1])
                tlen = int(row[6])
                qstart, qend = int(row[2]), int(row[3])
                tstart, tend = int(row[7]), int(row[8])
                strand = row[4]
                print('parsing', qstart, qend, tstart, tend)
                if strand == '+':
                    plt.plot([qstart, qend], [tstart, tend], color='red')
                else:
                    plt.plot([qstart, qend], [tend, tstart], color='blue')

    plt.xlabel('query: {}'.format(args.query))
    plt.ylabel('target: {}'.format(args.target))

    if args.qs or args.qe:
        plt.xlim(args.qs, args.qe)
    else:
        plt.xlim(0, qlen)

    if args.ts or args.te:
        plt.ylim(args.ts, args.te)
    else:
        plt.ylim(0, tlen)

    if args.out:
        plt.savefig(args.out)
    else:
        plt.show()


if __name__ == '__main__':
    main()
