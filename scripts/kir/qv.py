#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
import argparse
from pathlib import Path
from pafpy import PafFile


def main():
    parser = argparse.ArgumentParser(description='calculate QV')
    parser.add_argument('paf', type=Path, help='PAF')
    args = parser.parse_args()

    n = 0
    m = 0

    with PafFile(args.paf) as paf:
        for record in paf:
            if record.is_primary():
                qname = record.qname
                qlen = record.qlen
                n += qlen
                n_mismatch = record.get_tag("NM").value
                m += n_mismatch
                print('record', qname, qlen, n_mismatch)

    print('n_bases={} n_mismatch={} {}%'.format(n, m, m / n * 100))


if __name__ == '__main__':
    main()
