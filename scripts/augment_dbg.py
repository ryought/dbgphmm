#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import argparse
from pathlib import Path
import networkx as nx


def main():
    parser = argparse.ArgumentParser(description='Augment GFA of DBG')
    parser.add_argument('gfa', type=Path, help='GFA')
    parser.add_argument('k', type=int, help='k')
    args = parser.parse_args()
    k = args.k

    graph = nx.DiGraph()

    segments = []
    links = []

    with open(args.gfa) as f:
        for line in f:
            columns = line.split('\t')
            if columns[0] == 'S':
                id = columns[1]
                seq = columns[2]
                graph.add_node(id, seq=seq)
                segments.append(columns)
            elif columns[0] == 'L':
                id_a = columns[1]
                id_b = columns[3]
                graph.add_edge(id_a, id_b)
                links.append(columns)

    for segment in segments:
        id = segment[1]
        seq = segment[2]
        if graph.in_degree(id) > 0:
            prefixes = [graph.nodes[pred]['seq'][-(k-1):]
                        for pred in graph.predecessors(id)]
            assert all(prefix == prefixes[0] for prefix in prefixes)
            segment[2] = prefixes[0] + seq
        print('\t'.join(segment), end='')
    for link in links:
        print('\t'.join(link), end='')


if __name__ == '__main__':
    main()
