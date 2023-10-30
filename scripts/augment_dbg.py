#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import argparse
from pathlib import Path
import networkx as nx


def pred(graph, node, length):
    v = node
    s = ''
    while len(s) < length:
        v = next(graph.predecessors(v), None)
        if v is None:
            return s
        s = graph.nodes[v]['seq'] + s
    return s[-length:]


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
        prefix = pred(graph, id, k-1)
        segment[2] = prefix + seq
        # print(id, len(prefix), len(seq))
        print('\t'.join(segment), end='')
    for link in links:
        print('\t'.join(link), end='')


if __name__ == '__main__':
    main()
