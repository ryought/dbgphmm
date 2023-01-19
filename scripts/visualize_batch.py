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
import matplotlib.pyplot as plt
from dataclasses import dataclass
from visualize_posterior import parse_prob
import pickle

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
    return 'U{}N{}H{}P{}p{}'.format(U, N, H, str(P).replace('.', ''), str(p).replace('.', ''))

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('log_dir', type=Path, help='batch log directory')
    args = parser.parse_args()

    kmers = dict()
    samples = dict()
    for U in [1000, 500, 100, 50, 20]:
        N = 1000 // U
        P = 2
        # for H in [0.001, 0.005]:
        H = 0.001
        for p in [0.001, 0.005, 0.01]:
            kmers, samples, opts = parse_log_file(args.log_dir / filename_from_params(U, N, H, P, p))
            kmers[(U, H, p)] = kmers
            samples[(U, H, p)] = samples
    with open('batch_logs.pkl', 'wb') as f:
        pickle.dump((kmers, samples), f)

if __name__ == '__main__':
    main()
