#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from Bio import SeqIO
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description='Join fasta records')
    parser.add_argument('fasta', type=Path, help='FASTA')
    parser.add_argument('query', type=str,
                        help='record ids delimited by comma')
    parser.add_argument('--name', type=str,
                        default='out',
                        help='record ids delimited by comma')
    args = parser.parse_args()
    fasta = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))
    seq = ''
    for id in args.query.split(','):
        if id not in fasta:
            print('"{}" is not in fasta'.format(id))
            return
        seq += fasta[id].seq

    print('>{}'.format(args.name))
    print(seq)


if __name__ == '__main__':
    main()
