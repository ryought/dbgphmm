#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import argparse
import json
import csv
from dataclasses import dataclass
from enum import Enum
import numpy as np
import matplotlib
# matplotlib.use("TkAgg")
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import logsumexp

class PlotMode(Enum):
    Matrix = 1
    InitProbDistribution = 2
    ActiveNodeRatio = 3
    SortedDistribution = 4

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

    mode = PlotMode.SortedDistribution

    if mode == PlotMode.InitProbDistribution:
        for i in range(1, 10):
            plt.plot(probs[:, i])
            plt.savefig(args.tsv_filename + '.p{}.png'.format(i), dpi=200)
            plt.clf()
            plt.close()
    elif mode == PlotMode.Matrix:
        # for y in range(n_bases):
        #     print(y, np.exp(np.max(probs[:, y])))
        plt.figure(figsize=(30, 8))
        plt.subplot(2, 1, 1)
        plt.imshow(np.exp(probs.T), vmin=0, vmax=0.3)
        plt.colorbar()
        plt.subplot(2, 1, 2)
        plt.imshow(probs.T, vmin=-50, vmax=0)
        # plt.imshow(probs.T)
        plt.colorbar()
        plt.savefig(args.tsv_filename + '.png', dpi=200)
    elif mode == PlotMode.ActiveNodeRatio:
        # n_active_nodes_choice = [2, 8, 32, 128]
        n_active_nodes_choice = [2, 10, 50, 250]
        plt.figure(figsize=(30, 8))
        ratio = np.zeros((len(n_active_nodes_choice), n_bases))
        for t, n_active_nodes in enumerate(n_active_nodes_choice):
            for i in range(n_bases):
                ratio[t, i] = np.exp(logsumexp(np.sort(probs[:, i])[::-1][:n_active_nodes]))
                # print(ratio[i])
            plt.plot(ratio[t], label='n={}'.format(n_active_nodes))
        plt.grid(axis='y')
        plt.ylim(0, 1.1)
        plt.legend()
        plt.savefig(args.tsv_filename + '.activenoderatio.png', dpi=200)
    elif mode == PlotMode.SortedDistribution:
        plt.figure(figsize=(30, 8))
        is_log = False
        n_ignore_start = 16
        if is_log:
            plt.ylim(-100, 0)
            for i in range(n_ignore_start, n_bases):
                plt.plot(np.sort(probs[:, i])[::-1])
        else:
            plt.ylim(0 - 0.1, 1 + 0.1)
            for i in range(n_ignore_start, n_bases):
                plt.plot(np.exp(np.sort(probs[:, i])[::-1]))
        plt.savefig(args.tsv_filename + '.sorteddist.png', dpi=200)
    # plt.show()

if __name__ == '__main__':
    main()
