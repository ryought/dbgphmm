#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
python3 em_plotter.py hoge.tsv hoge.json
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import glob

from parsers import *

def kmer_existence_error(copy_nums, copy_nums_true):
    return sum([
        0 if ((y > 0 and x > 0) or (y == 0 and x == 0)) else 1
        for (x, y) in zip(copy_nums, copy_nums_true)
    ])

def missing_true_kmer_count(copy_nums, copy_nums_true):
    return sum([
        1 if (y > 0 and x == 0) else 0
        for (x, y) in zip(copy_nums, copy_nums_true)
    ])

def existing_false_kmer_count(copy_nums, copy_nums_true):
    return sum([
        1 if (y == 0 and x > 0) else 0
        for (x, y) in zip(copy_nums, copy_nums_true)
    ])

def diff(copy_nums, copy_nums_true):
    return sum([abs(x - y) for (x, y) in zip(copy_nums, copy_nums_true)])

def plot_history(logs, stats, filename):
    copy_nums_true = [log.copy_nums for log in logs if log.id == 'true'][0]

    N = 2

    plt.subplot(N, 1, 1)
    x = [missing_true_kmer_count(log.copy_nums, copy_nums_true) for log in logs if log.id != 'true']
    y = [existing_false_kmer_count(log.copy_nums, copy_nums_true) for log in logs if log.id != 'true']
    z = [diff(log.copy_nums, copy_nums_true) for log in logs if log.id != 'true']
    plt.plot(x, label='missing_true_kmer({})'.format(x[-1]))
    plt.plot(y, label='existing_false_kmer({})'.format(y[-1]))
    plt.plot(z, label='diff({})'.format(z[-1]))
    plt.legend()

    plt.subplot(N, 1, 2)
    plt.plot([log.depth for log in logs if log.id != 'true'])
    plt.ylabel('expected depth')

    plt.xlabel('iteration')
    plt.savefig(filename)

def plot_history_agg(logs_dict, stats, filename):
    for tag, logs in logs_dict.items():
        copy_nums_true = [log.copy_nums for log in logs if log.id == 'true'][0]
        x = [missing_true_kmer_count(log.copy_nums, copy_nums_true) for log in logs if log.id != 'true']
        y = [existing_false_kmer_count(log.copy_nums, copy_nums_true) for log in logs if log.id != 'true']
        z = [diff(log.copy_nums, copy_nums_true) for log in logs if log.id != 'true']
        plt.plot(x, label='missing_true_kmer({})'.format(x[-1]))
        plt.plot(y, label='existing_false_kmer({})'.format(y[-1]))
        plt.plot(z, label='diff({})'.format(z[-1]))
    plt.savefig(filename)

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('json_filename', type=str, help='json created by dbgphmm stats')
    parser.add_argument('tsv_filenames', nargs='+', type=str, help='')
    args = parser.parse_args()

    # # parse
    # stats = parse_stats(args.json_filename)
    # logs_dict = dict()
    # for tsv_filename in args.tsv_filenames:
    #     tag = parse_run_tag(tsv_filename)
    #     logs = parse_logs(tsv_filename)
    #     if logs:
    #         logs_dict[tag] = logs

    # # foreach
    # for tsv_filename in args.tsv_filenames:
    #     tag = parse_run_tag(tsv_filename)
    #     plot_history(logs_dict[tag], stats, tsv_filename + '.history.png')

    # collective
    plot_history_agg(logs_dict, stats, '')

def main_collective():
    for tsv_filename in glob.glob('./data/L100_N1_S*/E0.01_D10/fullemlingrad2_K*.tsv'):
        tag = parse_run_tag(tsv_filename)
        logs = parse_em_logs(tsv_filename)
        print(logs)

    for tsv_filename in glob.glob('./data/L100_N1_S*/E0.01_D10/freqem2_K*.tsv'):
        tag = parse_run_tag(tsv_filename)
        logs, _ = parse_logs(tsv_filename)
        print(logs)

if __name__ == '__main__':
    # main()
    main_collective()
