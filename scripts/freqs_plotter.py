#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
visualize kmer frequency history
python3 freqs_plotter.py hoge.txt hoge.json
"""
import argparse
import json
import csv
from dataclasses import dataclass
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from plotter import parse_stats

@dataclass
class FreqsLog:
    freqs: list

def parse_logs(txt_filename):
    with open(txt_filename) as f:
        logs = []
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            log = FreqsLog(
                freqs=json.loads(row[0][6:])
            )
            logs.append(log)
    return logs

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('json_filename', type=str, help='json created by dbgphmm stats')
    parser.add_argument('txt_filename', type=str, help='txt created by dbgphmm benchmark float-em')
    args = parser.parse_args()

    # parse
    stats = parse_stats(args.json_filename)
    logs = parse_logs(args.txt_filename)
    print(logs)

if __name__ == '__main__':
    main()
