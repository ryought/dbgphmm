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
from dataclasses import dataclass
from visualize_posterior import parse_prob

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

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('log_filename', type=str, help='batch log file')
    args = parser.parse_args()

    kmers, samples, opts = parse_log_file(args.log_filename)

if __name__ == '__main__':
    main()
