#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""
from pafpy import PafFile
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('paf', type=Path, help='PAF')
    args = parser.parse_args()

    with PafFile(args.paf) as paf:
        for record in paf:
            print(record.tname)
