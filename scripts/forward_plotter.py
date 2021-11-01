#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import argparse
import json
import csv
from dataclasses import dataclass
import numpy as np
import matplotlib
# matplotlib.use("TkAgg")
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import logsumexp

def parse_matrix(tsv_filename):
    kmers = []
    probs = []
    with open(tsv_filename) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            kmers.append(row[0])
            probs.append(list(map(lambda e: float(e), row[1:])))
    probs = np.array(probs)
    return kmers, probs

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('tsv_filename', type=str, help='tsv created by dbgphmm forward')
    args = parser.parse_args()

    # parse
    kmers, probs = parse_matrix(args.tsv_filename)
    n_kmer, n_bases = probs.shape

    for y in range(n_bases):
        probs[:, y] = probs[:, y] - logsumexp(probs[:, y])

    plt.figure(figsize=(30, 8))
    plt.subplot(2, 1, 1)
    plt.imshow(np.exp(probs.T))
    plt.subplot(2, 1, 2)
    plt.imshow(probs.T)
    # plt.imshow(probs.T)

    plt.savefig(args.tsv_filename + '.png', dpi=200)
    # plt.show()

if __name__ == '__main__':
    main()
