#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualize difference of MAP file (dbgphmm)
"""
import csv
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict
import argparse
from pathlib import Path
from dataclasses import dataclass
import re
import numpy as np
import glob
csv.field_size_limit(1000000000)


@dataclass(frozen=True)
class Emit:
    node: int
    prob: float


def parse_emits_list(s: str) -> List[Emit]:
    """
    >>> a = parse_emits_list('42:99.1,3:0.8,1:0.0,24:0.0')
    >>> a == [Emit(node=42, prob=99.1), Emit(node=3, prob=0.8), Emit(node=1, prob=0.0), Emit(node=24, prob=0.0)]
    True
    """
    return [Emit(node=int(t.split(':')[0]), prob=float(t.split(':')[1]))
            for t in s.split(',')]


def parse_map_file(filename: Path) -> List[List[List[Emit]]]:
    """
    """
    ret = []
    current_read = None
    currens_pos = None
    with open(filename) as f:
        for line in f:
            if line[0] == '#':
                continue
            s = line.rstrip().split('\t')
            if s:
                read = int(s[0])
                pos = int(s[1])
                base = s[2]
                emits = parse_emits_list(s[3])

                if current_read != read:
                    print(current_read)
                    ret.append([])
                ret[read].append(emits)
                current_read = read
                current_pos = pos
    return ret


def jaccard(a: List[Emit], b: List[Emit]) -> float:
    """
    >>> a = [Emit(node=42, prob=99.1), Emit(node=3, prob=0.8), Emit(node=1, prob=0.0)]
    >>> b = [Emit(node=42, prob=99.1), Emit(node=3, prob=0.8), Emit(node=0, prob=0.0)]
    >>> jaccard(a, b) == 2/4
    True
    """
    sa = set(e.node for e in a)
    sb = set(e.node for e in b)
    return len(sa & sb) / len(sa | sb)


def jaccards(map_a, map_b):
    ret = []
    n_reads = len(map_a)
    for i in range(n_reads):
        s = [jaccard(map_a[i][j], map_b[i][j]) for j in range(len(map_a[i]))]
        ret.append(s)
    return ret


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('map_a', type=Path, help='MAP file')
    parser.add_argument('map_b', type=Path, help='MAP file')
    # parser.add_argument('--query', action='store_true')
    args = parser.parse_args()

    map_a = parse_map_file(args.map_a)
    map_b = parse_map_file(args.map_b)
    assert(len(map_a) == len(map_b))
    n = len(map_a)
    for i in range(n):
        assert(len(map_a[i]) == len(map_b[i]))

    js = jaccards(map_a, map_b)
    for j in range(len(js[0])):
        print(j, js[0][j])
    # for i in range(10):
    #     plt.plot(js[i])
    # plt.show()


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    main()
