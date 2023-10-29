#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from Bio import SeqIO
import argparse
from pathlib import Path
import sys
import json
from collections import defaultdict

n_samples = 0
total_time = 0

d = Path('sim/n4_p0.0003/H0.01_H00.0001/dbgphmm')
for post in d.glob('*.post'):
    print(post)
    n_samples_of_k = 0
    total_time_of_k = 0
    with post.open() as f:
        for line in f:
            if line[0] == 'C':
                body = line.split('\t')[3]
                t = json.loads(body)["time_likelihood"]
                n_samples += 1
                total_time += t
                n_samples_of_k += 1
                total_time_of_k += t
                # print(t)
    print(post.name, n_samples_of_k, total_time_of_k, total_time_of_k / 32)
print('n_samples', n_samples)
print('total_time', total_time)
print('average_time', total_time / n_samples)
