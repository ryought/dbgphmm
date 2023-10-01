#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from Bio import SeqIO
import argparse
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser(description='Join fasta records')
    parser.add_argument('fasta', type=Path, help='FASTA')
    parser.add_argument('query', type=str, nargs='+',
                        help='record ids delimited by comma and space')
    parser.add_argument('--name', type=str,
                        default='out',
                        help='')
    args = parser.parse_args()
    fasta = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

    for i, query in enumerate(args.query):
        seq = ''
        for q in query.split(','):
            direction = q[-1]
            id = q[:-1]
            if id not in fasta:
                print('"{}" is not in fasta'.format(id), file=sys.stderr)
                sys.exit(1)
            if direction == '+':
                seq += fasta[id].seq
            elif direction == '-':
                seq += fasta[id].seq.reverse_complement()
            else:
                print('invalid direction "{}"'.format(
                    direction), file=sys.stderr)
                sys.exit(1)
        print('>{}{}'.format(args.name, i))
        print(seq)


if __name__ == '__main__':
    main()
