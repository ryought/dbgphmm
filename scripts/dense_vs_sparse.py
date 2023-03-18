#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import csv
csv.field_size_limit(1000000000)
from collections import defaultdict
from parse import *
import math
import statistics
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('filename', type=str, help='')
    args = parser.parse_args()

    # parse
    opts = ''
    reads = defaultdict(list)

    with open(args.filename) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row[0][0] == '#':
                continue
            read_id = int(row[0])
            # read_pos = int(row[1])
            p = float(row[4])
            reads[read_id].append(p)

    for (read_id, read) in reads.items():
        plt.plot(read, label='{}'.format(read_id))
    plt.legend()
    plt.show()

main()
