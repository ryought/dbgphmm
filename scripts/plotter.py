#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
visualize training log

python3 plotter.py hoge.csv hoge.json
"""
import argparse
import json
import csv
from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt

@dataclass
class OptimizeLog:
    id: int
    temp: float
    now_score: float
    next_score: float
    p_accept: float
    is_accepted: bool
    now_score_prior: float
    now_score_forward: float
    now_size: int
    now_state: list
    next_score_prior: float
    next_score_forward: float
    next_size: int
    next_state: list

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('tsv_filename', type=str, help='tsv created by dbgphmm optimize')
    parser.add_argument('json_filename', type=str, help='json created by dbgphmm stats')
    args = parser.parse_args()

    with open(args.tsv_filename) as f:
        logs = []
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            log = OptimizeLog(
                id=int(row[0]),
                temp=float(row[1]),
                now_score=float(row[2]),
                next_score=float(row[3]),
                p_accept=float(row[4]),
                is_accepted=(row[5] == 'true'),
                now_score_prior=float(row[6]),
                now_score_forward=float(row[7]),
                now_size=int(row[8]),
                now_state=json.loads(row[9]),
                next_score_prior=float(row[10]),
                next_score_forward=float(row[11]),
                next_size=int(row[12]),
                next_state=json.loads(row[13]),
            )
            logs.append(log)

    with open(args.json_filename) as f:
        stats = json.load(f)

if __name__ == '__main__':
    main()
