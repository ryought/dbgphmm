#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import argparse
import csv
from collections import defaultdict
from parse import *
from pathlib import Path
import math
import numpy as np
from itertools import product
import matplotlib.pyplot as plt
from dataclasses import dataclass
from visualize_posterior import parse_prob
import pickle

Us = [1000, 500, 100, 50, 20]
L = 1000
P = 2
Hs = [0.001, 0.005]
H = 0.001
ps = [0.001, 0.005, 0.01]
markers = ['+', '*', 'x']

@dataclass
class Experiment:
    k: int
    seed: int

@dataclass
class Sample:
    post: float
    log_likelihood: float
    g: int
    diff: int
    n_missing: int
    info: str

@dataclass
class Kmer:
    kmer: str
    node: int
    copy_num: int
    copy_num_true: int
    p_true: float
    p_zero: float

def parse_log_file(filename):
    kmers = defaultdict(lambda: defaultdict(list))
    samples = defaultdict(lambda: defaultdict(list))
    opts = ''
    with open(filename) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0].startswith('K'):
                k = int(row[1])
                seed = int(row[12])
                kmer = Kmer(
                    kmer=row[2],
                    node=int(row[3]),
                    copy_num=int(row[4]),
                    copy_num_true=int(row[5]),
                    p_true=float(row[7]),
                    p_zero=float(row[8]),
                )
                kmers[seed][k].append(kmer)
            elif row[0].startswith('N'):
                k = int(row[1])
                seed = int(row[2])
                sample = Sample(
                    post=parse_prob(row[3])[1],
                    log_likelihood=float(row[4]),
                    g=int(row[6]),
                    diff=int(row[9]),
                    n_missing=int(row[11]),
                    info=row[13],
                )
                samples[seed][k].append(sample)
            elif row[0].startswith('# opts='):
                opts = row[0]
            elif row[0].startswith('# '):
                pass
            else:
                pass
    return kmers, samples, opts

def filename_from_params(U, N, H, P, p):
    return 'U{}N{}H{}P{}p{}'.format(U, N, str(H).replace('.', ''), P, str(p).replace('.', ''))

def defaultdict_to_dict(d):
    for k, v in d.items():
        if isinstance(v, dict):
            d[k] = defaultdict_to_dict(v)
    return dict(d)

def fig_1(kmers, samples, filename=None):
    """
    kmer mis-classification as zero-copy rate

    x-axis: k
    y-axis: P(c[v]=0) s.t. c_true[v]=1
    """
    plt.figure(figsize=(10, 10))
    for i, (U, p) in enumerate(product(Us, ps)):
        plt.subplot(len(Us), len(ps), i + 1)
        for (seed, d) in kmers[(U, H, p)].items():
            xs = []
            ys = []
            for (k, dd) in d.items():
                for kmer in dd:
                    if kmer.copy_num_true == 1:
                        xs.append(k)
                        ys.append(kmer.p_zero)
            plt.scatter(xs, ys, marker=markers[seed], alpha=0.5)
        plt.ylim(0, 1)
        plt.xlim(10, 50)
        plt.xlabel('k')
        plt.ylabel('P(c[v]=0) c*[v]=1')
        plt.title('U={} p={}'.format(U, p))
    if filename:
        plt.savefig(filename, dpi=200)
    else:
        plt.show()


def fig_2(kmers, samples, filename=None):
    """
    missing k-mer count and model likelihood

    x-axis: k
    y-axis: P(c[v]=0) s.t. c_true[v]=1
    """
    # Us = [500]
    # ps = [0.01]
    plt.figure(figsize=(10, 10))
    plt.suptitle('Log likelihood difference between best models w/ and w/o missing kmers')
    for i, (U, p) in enumerate(product(Us, ps)):
        plt.subplot(len(Us), len(ps), i + 1)
        plt.title('U={} p={}'.format(U, p))
        plt.ylim(-10, 10)
        plt.grid()
        for seed in samples[(U, H, p)].keys():
            xs = []
            ys = []
            for (k, dd) in samples[(U, H, p)][seed].items():
                x0 = max(sample.log_likelihood for sample in dd if sample.n_missing == 0)
                x1 = max(sample.log_likelihood for sample in dd if sample.n_missing > 0)
                xs.append(k)
                ys.append(x0 - x1)
            plt.scatter(xs, ys)
    if filename:
        plt.savefig(filename, dpi=200)
    else:
        plt.show()
    plt.clf()
    plt.close()


def fig_supp2(kmers, samples, filename=None):
    """
    missing k-mer count and model likelihood
    show the distribution of missing-kmer and likelihood of sampled models
    for each seed and k value.

    x-axis: n_missings
    y-axis: log_likelihood
    """
    # Us = [500]
    # ps = [0.01]
    for _, (U, p) in enumerate(product(Us, ps)):
        plt.figure(figsize=(20, 10))
        plt.suptitle('n_missing vs log_likelihood (U={} p={})'.format(U, p))
        seeds = list(samples[(U, H, p)].keys())
        ks = list(samples[(U, H, p)][seeds[0]].keys())
        for i, (seed, k) in enumerate(product(seeds, ks)):
                plt.subplot(len(seeds), len(ks), i + 1)
                xs = []
                ys = []
                for sample in samples[(U, H, p)][seed][k]:
                    xs.append(sample.n_missing)
                    ys.append(sample.log_likelihood)
                plt.scatter(xs, ys, alpha=0.5)
                plt.xlim(0, 100)
                plt.title('seed={} k={}'.format(seed, k))

        if filename:
            plt.savefig(filename.format(U=U, p=p))
        else:
            plt.show()
        plt.clf()
        plt.close()

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('log_dir', type=Path, help='batch log directory')
    parser.add_argument('--to-pickle', type=Path, help='batch log directory')
    args = parser.parse_args()

    kmers = dict()
    samples = dict()
    for U in Us:
        N = L // U
        for p in ps:
            k, s, o = parse_log_file(args.log_dir / filename_from_params(U, N, H, P, p))
            kmers[(U, H, p)] = defaultdict_to_dict(k)
            samples[(U, H, p)] = defaultdict_to_dict(s)

    if args.to_pickle:
        with open(args.to_pickle, 'wb') as f:
            pickle.dump((kmers, samples), f)

    fig_1(kmers, samples, filename='fig_1.png')
    fig_2(kmers, samples, filename='fig_2.png')
    fig_supp2(kmers, samples, filename='fig_supp2_U{U}p{p}.png')

if __name__ == '__main__':
    main()
