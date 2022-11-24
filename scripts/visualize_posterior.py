#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import argparse
import csv
from collections import defaultdict
from parse import *
import math
import numpy as np
import matplotlib.pyplot as plt

def parse_post(s):
    """
    parse posterior distribution of kmer copy numbers
    """
    ret = defaultdict(float)
    for p in s.split(','):
        x, log_p_x, p_x = parse('p(x={})={}({})', p)
        x = int(x)
        log_p_x = float(log_p_x)
        ret[x] = math.exp(log_p_x)
    return ret

def test_parse_post():
    print(parse_post('p(x=0)=-2683.066794904139(0.0000),p(x=1)=-0.0011576044833668843(0.9988),p(x=2)=-6.761981256076189(0.0012)'))

def parse_prob(s):
    log_p_x, p_x = parse('{}({})', s)
    log_p_x = float(log_p_x)
    p_x = math.exp(log_p_x)
    return (log_p_x, p_x)

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('log_filename', type=str, help='sample_posterior log file')
    args = parser.parse_args()

    # parse
    opts = ''
    n_nodes_of_k = dict()
    n_edges_of_k = dict()
    n_purged_of_k = dict()
    kmer_post_of_k = defaultdict(list)
    post_of_k = defaultdict(list)
    with open(args.log_filename) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0].startswith('K'):
                k = int(row[1])
                copy_num_true = int(row[4])
                post = parse_post(row[5])
                kmer_post_of_k[k].append((copy_num_true, post[copy_num_true]))
            elif row[0].startswith('N'):
                k = int(row[1])
                lp, p = parse_prob(row[2])
                ll, l = parse_prob(row[3])
                g = int(row[5])
                d = int(row[7])
                post_of_k[k].append((p, g, d, ll))
                pass
            elif row[0].startswith('# opts='):
                opts = row[0]
            elif row[0].startswith('# '):
                seg = row[0].split(' ')
                if len(seg) == 3:
                    k = int(seg[1].lstrip('k='))
                    # n_nodes
                    if seg[2].startswith('n_nodes='):
                        n_nodes_of_k[k] = int(seg[2].lstrip('n_nodes='))
                    # n_edges
                    if seg[2].startswith('n_edges='):
                        n_edges_of_k[k] = int(seg[2].lstrip('n_edges='))
                    # n_purged
                    if seg[2].startswith('n_purged='):
                        n_purged_of_k[k] = int(seg[2].lstrip('n_purged='))
                pass
            else:
                pass
    print(opts)
    print(n_nodes_of_k)
    print(n_edges_of_k)
    print(n_purged_of_k)

    # visualize
    # (1) graph structure
    plt.figure(figsize=(20, 10))
    plt.subplot(4, 1, 1)
    ks = list(n_nodes_of_k.keys())
    n_nodes = list(n_nodes_of_k.values())
    plt.plot(ks, n_nodes, label='n_nodes', marker='o')
    ks = list(n_edges_of_k.keys())
    n_edges = list(n_edges_of_k.values())
    plt.plot(ks, n_edges, label='n_edges', marker='o')
    ks = list(n_purged_of_k.keys())
    n_purged = list(n_purged_of_k.values())
    plt.plot(ks, n_purged, label='n_purged', marker='o')
    plt.xticks(ks, fontsize=8)
    plt.ylabel('graph props')
    plt.legend()

    # (2) posterior distribution
    plt.subplot(4, 1, 2)
    ks_true = []
    ps_true = []
    for k, post in post_of_k.items():
        ps = [p for (p, g, d, ll) in post if d != 0]
        ks = [k for (p, g, d, ll) in post if d != 0]
        plt.scatter(ks, ps, c='blue', marker='x', alpha=0.5)

        for (p, g, d, ll) in post:
            if d == 0:
                ks_true.append(k)
                ps_true.append(p)

    print(len(ks_true), len(ps_true))
    plt.plot(ks_true, ps_true, c='red', marker='o', alpha=0.5)
    plt.xticks(ks_true, fontsize=8)
    plt.ylabel('P(G|R)')

    # (3) log likelihood
    plt.subplot(4, 1, 3)
    ks_true = []
    lls_true = []
    for k, post in post_of_k.items():
        lls = [ll for (p, g, d, ll) in post if d != 0]
        ks = [k for (p, g, d, ll) in post if d != 0]
        plt.scatter(ks, lls, c='blue', marker='x', alpha=0.5)
        for (p, g, d, ll) in post:
            if d == 0:
                ks_true.append(k)
                lls_true.append(ll)
    plt.plot(ks_true, lls_true, c='red', marker='o', alpha=0.5)
    plt.xticks(ks_true, fontsize=8)
    plt.ylabel('log P(R|G)')

    # (4) copy_num posterior
    plt.subplot(4, 1, 4)
    plt.ylabel('P(c[v] = c_true[v])')
    plt.xlabel('k')

    plt.suptitle(opts, wrap=True)
    plt.show()

if __name__ == '__main__':
    main()
    # test_parse_post()
