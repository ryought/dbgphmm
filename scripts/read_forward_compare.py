#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import argparse
from dataclasses import dataclass
from enum import Enum
from parse import *
import numpy as np
import matplotlib.pyplot as plt
import csv
from collections import defaultdict

@dataclass
class ReadForward:
    read_id: int
    base: str
    p_a: float
    p_b: float
    dp_a: float
    dp_b: float
    ddp: float

def plot_history(rs):
    n = len(rs.keys())
    plt.figure(figsize=(5, 20))
    for (i, read_id) in enumerate(rs.keys()):
        plt.subplot(n, 1, i + 1)
        x = [r.base for r in rs[read_id]]
        p_a = [r.p_a for r in rs[read_id]]
        p_b = [r.p_b for r in rs[read_id]]
        plt.title('{}'.format(read_id))
        plt.plot(x, p_a)
        plt.plot(x, p_b)
    plt.savefig('hoge.pdf')

def plot_history_diff(rs, ro):
    plt.subplot(2, 1, 1)
    for (i, read_id) in enumerate(rs.keys()):
        if ro[read_id] == 0:
            x = [r.base for r in rs[read_id]]
            ddp = [r.ddp for r in rs[read_id]]
            plt.plot(x, ddp, label='{}'.format(read_id))
    plt.xticks(np.arange(0, 1100, 50))
    plt.legend()

    plt.subplot(2, 1, 2)
    for (i, read_id) in enumerate(rs.keys()):
        if ro[read_id] == 1:
            x = [r.base for r in rs[read_id]]
            ddp = [r.ddp for r in rs[read_id]]
            plt.plot(x, ddp, label='{}'.format(read_id))
    plt.xticks(np.arange(0, 1100, 50))
    plt.legend()

    plt.show()

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('summary_filename', type=str)
    args = parser.parse_args()

    ro = dict()
    rs = defaultdict(list)
    with open(args.summary_filename) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) == 0:
                pass
            else:
                if row[0].startswith('r['):
                    read_id = int(row[0][2:-1])
                    if row[1] == 'summary':
                        origin = int(row[2])
                        ro[read_id] = origin
                    else:
                        base = row[1]
                        p_a = float(row[3])
                        p_b = float(row[4])
                        dp_a = float(row[5])
                        dp_b = float(row[6])
                        ddp = float(row[7])
                        r = ReadForward(
                            read_id=read_id,
                            base=base,
                            p_a=p_a,
                            p_b=p_b,
                            dp_a=dp_a,
                            dp_b=dp_b,
                            ddp=ddp,
                        )
                        rs[read_id].append(r)

    # plot_history(rs)
    plot_history_diff(rs, ro)

if __name__ == '__main__':
    main()
