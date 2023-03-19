#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
csv.field_size_limit(1000000000)
from collections import defaultdict
from parse import *
import math
import statistics
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('filename', type=str, help='')
    parser.add_argument('--backward', action='store_true')
    parser.add_argument('--show', action='store_true')
    args = parser.parse_args()

    # parse
    opts = ''
    ps = defaultdict(list)
    nms = defaultdict(list)
    nis = defaultdict(list)
    nds = defaultdict(list)

    with open(args.filename) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0][0] == '#':
                continue
            read_id = int(row[0])
            read_pos = int(row[1])
            p = float(row[4])
            nm = int(row[7])
            ni = int(row[8])
            nd = int(row[9])
            ps[read_id].append(p)
            nms[read_id].append(nm)
            nis[read_id].append(ni)
            nds[read_id].append(nd)
            print(read_id, read_pos, p, nm, ni, nd)

    # for (read_id, read) in ps.items():
    #     plt.plot(read, label='{}'.format(read_id))
    plt.figure(figsize=(20, 10))
    plt.subplot(3, 1, 1)
    for read_id in nms.keys():
        if args.backward:
            nms[read_id].reverse()
        plt.plot(nms[read_id], label='m{}'.format(read_id))
    plt.legend()
    plt.subplot(3, 1, 2)
    for read_id in nis.keys():
        if args.backward:
            nis[read_id].reverse()
        plt.plot(nis[read_id], label='i{}'.format(read_id))
    plt.legend()
    plt.subplot(3, 1, 3)
    for read_id in nds.keys():
        if args.backward:
            nds[read_id].reverse()
        plt.plot(nds[read_id], label='d{}'.format(read_id))
    plt.legend()

    if args.show:
        plt.show()
    else:
        plt.savefig(args.filename + '.png', dpi=200)

main()
