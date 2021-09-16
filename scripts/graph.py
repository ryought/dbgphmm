#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import argparse
import networkx as nx
import csv
import matplotlib.pyplot as plt
from fa2 import ForceAtlas2

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('layout_filename', type=str, help='layout created by dbgphmm visualize')
    args = parser.parse_args()

    nodes, edges = dict(), dict()
    pos = dict()
    with open(args.layout_filename, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for i, row in enumerate(reader):
            if row[0] == 'N':
                # node
                km1mer, x, y = row[1], float(row[2]), float(row[3])
                nodes[km1mer] = (x, y)
            elif row[0] == 'E':
                # edge
                kmer, c = row[1], float(row[2])
                edges[kmer] = c

    G = nx.MultiDiGraph()
    for node in nodes.keys():
        G.add_node(node)
    for edge, count in edges.items():
        prefix = edge[:-1]
        suffix = edge[1:]
        G.add_edge(prefix, suffix)

    forceatlas2 = ForceAtlas2(verbose=True)
    pos = forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=10000)
    nx.draw(G, pos=pos)
    # nx.draw(G)

    plt.savefig(args.layout_filename + '.png', dpi=200)

if __name__ == '__main__':
    main()
