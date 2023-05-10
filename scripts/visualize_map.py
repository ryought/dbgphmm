#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Visualize difference of MAP file (dbgphmm)
"""
from dbgphmm import parse_dbg_file, parse_map_file, jaccards
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


def edge_map(dbg):
    ret = {}
    for (s, t, w) in dbg.edges(data=True):
        E = w['id']
        for (i, e) in enumerate(w['full_edges']):
            ret[e] = (E, i)
    return ret


def to_compact_edge(emits, edge_map):
    """
    hogefuga
    """
    return [edge_map[e.node] for e in emits]


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--map_a', type=Path, help='MAP file')
    parser.add_argument('--map_b', type=Path, help='MAP file')
    parser.add_argument('--dbg', type=Path, help='DBG file')
    # parser.add_argument('--query', action='store_true')
    args = parser.parse_args()

    dbg = parse_dbg_file(args.dbg)
    m = edge_map(dbg)
    print(m)

    if True:
        map_a = parse_map_file(args.map_a)
        map_b = parse_map_file(args.map_b)
        assert(len(map_a) == len(map_b))
        n = len(map_a)
        for i in range(n):
            assert(len(map_a[i]) == len(map_b[i]))

        js = jaccards(map_a, map_b)
        i = 0
        for j in range(len(js[i])):
            print(j, js[i][j])
            # sort by position in compacted graph
            # print(sorted(to_compact_edge(map_a[i][j], m)))
            # print(sorted(to_compact_edge(map_b[i][j], m)))
            # sort by estimated state probability
            print(to_compact_edge(map_a[i][j], m))
            print(to_compact_edge(map_b[i][j], m))
        # for i in range(10):
        #     plt.plot(js[i])
        # plt.show()


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    main()
