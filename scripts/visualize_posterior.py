#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import argparse
import csv
csv.field_size_limit(1000000000)
from collections import defaultdict
from parse import *
import math
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

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

@dataclass
class Sample:
    # posterior probability P(G|R)
    p: float
    # log posterior probability logP(G|R)
    lp: float
    # log likelihood logP(R|G)
    ll: float
    # genome size
    g: int
    # iteration id
    i: int
    # dist from true copy nums
    d: int

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('log_filename', type=str, help='sample_posterior log file')
    args = parser.parse_args()

    # parse
    opts = ''
    n_nodes_of_k = dict()
    n_edges_of_k = dict()
    n_purged_of_k = dict()
    kmer_post_of_k = defaultdict(lambda: defaultdict(list))
    non_zero_kmer_post_of_k = defaultdict(lambda: defaultdict(list))
    post_of_k = defaultdict(list)
    with open(args.log_filename) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0] == 'K':
                k = int(row[1])
                copy_num_true = int(row[5])
                post = parse_post(row[9])
                kmer_post_of_k[copy_num_true][k].append(post[copy_num_true])
                if copy_num_true > 0:
                    non_zero_kmer_post_of_k[copy_num_true][k].append(post[0])
            elif row[0] == 'N':
                k = int(row[1])
                lp, p = parse_prob(row[2])
                ll = float(row[3])
                g = int(row[5])
                i = int(row[6])  # iteration id
                d = int(row[7])  # greedy search depth
                post_of_k[k].append(Sample(
                    lp=lp,
                    p=p,
                    ll=ll,
                    g=g,
                    i=i,
                    d=d,
                ))
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
    h = 3
    w = 3
    fig, ax = plt.subplots(h, w, sharex="all", figsize=(20, 10))
    ks_for_tick = list(n_nodes_of_k.keys())
    plt.xticks(ks_for_tick, fontsize=8)
    # (1) graph structure
    def draw_graph_structure(ax):
        ks = list(n_nodes_of_k.keys())
        n_nodes = list(n_nodes_of_k.values())
        ax.plot(ks, n_nodes, label='n_nodes', marker='o')
        ks = list(n_edges_of_k.keys())
        n_edges = list(n_edges_of_k.values())
        ax.plot(ks, n_edges, label='n_edges', marker='o')
        ks = list(n_purged_of_k.keys())
        n_purged = list(n_purged_of_k.values())
        ax.plot(ks, n_purged, label='n_purged', marker='x', color='black')
        twin1 = ax.twinx()
        twin1.bar(ks, n_purged, color='black', alpha=0.5)
        twin1.set_ylim(0, 10)
        twin1.set_ylabel('n_purged')
        ax.set_ylabel('graph props')
        ax.set_ylim(0, None)
        ax.legend()
    draw_graph_structure(ax[0, 0])

    def set_yaxis_as_prob_distribution(ax):
        ax.set_ylim(0-0.1, 1+0.1)
        ax.grid(axis='both')

    # (2) posterior distribution
    def draw_posterior_distribution(ax):
        ks_true = []
        ps_true = []
        for k, post in post_of_k.items():
            ps = [s.p for s in post if s.d != 0]
            ks = [k for s in post if s.d != 0]
            ax.scatter(ks, ps, c='blue', marker='x', alpha=0.5)

            for s in post:
                if s.d == 0:
                    ks_true.append(k)
                    ps_true.append(s.p)

        print(len(ks_true), len(ps_true))
        ax.plot(ks_true, ps_true, c='red', marker='o', alpha=0.5)
        set_yaxis_as_prob_distribution(ax)
        ax.set_ylabel('P(G|R)')
    draw_posterior_distribution(ax[1, 0])

    # (3) log likelihood
    def draw_log_likelihood(ax):
        ks = []
        lls = []
        gs = []
        ks_true = []
        lls_true = []
        for k, post in post_of_k.items():
            lls = lls + [s.ll for s in post if s.d != 0]
            ks = ks + [k for s in post if s.d != 0]
            gs = gs + [s.g for s in post if s.d != 0]
            # ax.scatter(ks, lls, s=gs, c='blue', marker='o', alpha=0.5)
            for s in post:
                if s.d == 0:
                    ks_true.append(k)
                    lls_true.append(s.ll)
        pcm = ax.scatter(ks, lls, c=gs, marker='o', alpha=0.3)
        fig.colorbar(pcm, ax=ax, location='bottom', label='genome_size')
        ax.plot(ks_true, lls_true, c='red', marker='*', alpha=0.5)
        ax.set_ylabel('log P(R|G)')
    draw_log_likelihood(ax[2, 0])

    # (4) copy_num posterior
    def draw_copy_num_posterior(ax):
        for copy_num_true, post_of_k in kmer_post_of_k.items():
            ks = [k for k, kmer_post in post_of_k.items() for _ in kmer_post]
            ps = [p for k, kmer_post in post_of_k.items() for p in kmer_post]
            ax.scatter(ks, ps, marker='o', alpha=0.5, label='x{}'.format(copy_num_true))
        ax.legend()
        ax.set_ylabel('P(c[v]=c_true[v])')
        set_yaxis_as_prob_distribution(ax)
    draw_copy_num_posterior(ax[0, 1])

    # (5) copy_num posterior of 0x nodes
    def draw_copy_num_posterior_0x(ax):
        copy_num_true = 0
        post_of_k = kmer_post_of_k[copy_num_true]
        ks = [k for k, kmer_post in post_of_k.items() for _ in kmer_post]
        ps = [p for k, kmer_post in post_of_k.items() for p in kmer_post]
        ax.scatter(ks, ps, marker='o', alpha=0.5, label='x{}'.format(copy_num_true))
        ax.legend()
        ax.set_ylabel('P(c[v]=c_true[v])')
        set_yaxis_as_prob_distribution(ax)
    draw_copy_num_posterior_0x(ax[1, 1])

    # (5) copy_num posterior of >0x nodes
    def draw_copy_num_posterior_non_0x(ax):
        for copy_num_true, post_of_k in non_zero_kmer_post_of_k.items():
            assert(copy_num_true != 0)
            ks = [k for k, kmer_post in post_of_k.items() for _ in kmer_post]
            ps = [p for k, kmer_post in post_of_k.items() for p in kmer_post]
            ax.scatter(ks, ps, marker='o', alpha=0.5, label='x{}'.format(copy_num_true))
        ax.legend()
        ax.set_ylabel('P(c[v]=0) v: c_true[v]!=0')
        set_yaxis_as_prob_distribution(ax)
    draw_copy_num_posterior_non_0x(ax[2, 1])

    # (6) n_sampled_copy_nums
    def draw_n_samples(ax):
        ks = [k for k, kmer_post in post_of_k.items()]
        ns = [len(kmer_post) for k, kmer_post in post_of_k.items()]
        ax.scatter(ks, ns, label='n_samples')
        ax.set_ylabel('n_samples')
        ax.set_ylim(0, None)
    draw_n_samples(ax[0, 2])

    # (6) n_sampled_copy_nums
    def draw_search_history(ax):
        ks = []
        lps = []
        Is = []
        for k, post in post_of_k.items():
            lps = lps + [s.lp for s in post]
            ks = ks + [k for s in post]
            Is = Is + [s.i for s in post]
        pcm = ax.scatter(ks, lps, c=Is, marker='o', alpha=0.2)
        fig.colorbar(pcm, ax=ax, location='bottom', label='iteration id')
        ax.set_ylabel('log P(G|R)')
        ax.grid(axis='x')
    draw_search_history(ax[1, 2])

    plt.suptitle(opts, wrap=True)
    # plt.show()
    plt.savefig(args.log_filename + '.summary.pdf')

if __name__ == '__main__':
    main()
    # test_parse_post()
