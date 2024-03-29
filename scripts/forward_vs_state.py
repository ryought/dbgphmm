#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dataclasses import dataclass
import matplotlib.pyplot as plt
import numpy as np
import statistics
import math
from parse import *
from collections import defaultdict
import argparse
import csv
csv.field_size_limit(1000000000)


def parse(s):
    """
    parse summary string of states vector

    >>> parse('N(v102):0.001,I(v1):0.9')
    [('N(v102)', 0.001), ('I(v1)', 0.9)]
    """
    return [(x.split(':')[0], float(x.split(':')[1])) for x in s.split(',')]


def jaccard(xs, ys):
    """
    calculate jaccard index

    >>> jaccard([('a'), ('b')], [('a'), ('c'), ('d')])
    0.25
    """
    sxs = set([x[0] for x in xs])
    sys = set([y[0] for y in ys])
    u = sxs | sys
    i = sxs & sys
    return len(i) / len(u)


def overlap(xs, ys):
    """
    how xs overs ys? Ratio of ys elements covered in xs

    >>> overlap([('a'), ('b')], [('a'), ('c')])
    0.5
    """
    sxs = set([x[0] for x in xs])
    sys = set([y[0] for y in ys])
    i = sxs & sys
    return len(i) / len(sys)


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('filename', type=str, help='')
    parser.add_argument('--show', action='store_true')
    parser.add_argument('--read-ids', nargs='*', type=int, help='')
    args = parser.parse_args()

    # parse
    forwards = defaultdict(list)
    states = defaultdict(list)
    js = defaultdict(list)

    with open(args.filename) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0][0] == '#':
                continue
            read_id = int(row[0])
            read_pos = int(row[1])
            forward = parse(row[4])
            state = parse(row[5])
            forwards[read_id].append(forward)
            states[read_id].append(state)
            # j = jaccard(forward, state)
            j = overlap(forward, state)
            js[read_id].append(j)

    # for (read_id, read) in ps.items():
    #     plt.plot(read, label='{}'.format(read_id))
    plt.figure(figsize=(20, 10))
    if args.read_ids:
        for read_id in args.read_ids:
            plt.plot(js[read_id], label='{}'.format(read_id))
    else:
        for (read_id, j) in js.items():
            plt.plot(j, label='{}'.format(read_id))
    plt.legend()
    plt.ylim(0, 1)

    if args.show:
        plt.show()
    else:
        plt.savefig(args.filename + '.png', dpi=200)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    main()
