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
    args = parser.parse_args()

    with open(args.paf) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0] == args.query and row[5] == args.target:
                qstart, qend = int(row[2]), int(row[3])
                tstart, tend = int(row[7]), int(row[8])
                print('parsing', qstart, qend, tstart, tend)
                plt.plot([qstart, qend], [tstart, tend])
    plt.xlabel('query')
    plt.ylabel('target')
    plt.show()


if __name__ == '__main__':
    main()
