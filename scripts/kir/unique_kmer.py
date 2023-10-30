#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""

from Bio import SeqIO
import argparse
from pathlib import Path
import sys
from collections import defaultdict


def disjoint_kmers(ref, asm, verbose=False):
    n_ref_only = 0
    n_asm_only = 0
    for kmer, occ in ref.items():
        if kmer not in asm:
            n_ref_only += 1
            if verbose:
                print('ref_only', kmer, occ)
    for kmer, occ in asm.items():
        if kmer not in ref:
            n_asm_only += 1
            if verbose:
                print('asm_only', kmer, occ)
    return n_ref_only, n_asm_only


def unique_kmers_in_ref(ref, asm):
    n_unique_kmers = 0
    n_unique_kmers_contained = 0
    for kmer, occ in ref.items():
        if len(occ) == 1:
            n_unique_kmers += 1
            if kmer in asm:
                n_unique_kmers_contained += 1
    return n_unique_kmers, n_unique_kmers_contained


def same_count(ref, asm, verbose=False):
    n_same = 0
    n_different = 0
    for kmer in (set(ref.keys()) | set(asm.keys())):
        if len(ref[kmer]) == len(asm[kmer]):
            n_same += 1
        else:
            n_different += 1
            if verbose:
                print('different', kmer, ref[kmer], asm[kmer])
    return n_same, n_different


def count_kmers_with_position(seqs, k):
    kmers = defaultdict(list)
    for h, seq in enumerate(seqs.values()):
        # forward
        s = sanitize(str(seq.seq))
        for i in range(len(s) - k + 1):
            kmer = s[i:i+k]
            kmers[kmer].append((h, i))
        # backward
        s = sanitize(str(seq.seq.reverse_complement()))
        for i in range(len(s) - k + 1):
            kmer = s[i:i+k]
            kmers[kmer].append((h, i))
    return kmers


def sanitize(seq: str):
    return seq.upper().replace('N', '')


def main():
    parser = argparse.ArgumentParser(description='unique k-mer analysis')
    parser.add_argument('ref', type=Path, help='FASTA')
    parser.add_argument('asm', type=Path, help='FASTA')
    parser.add_argument('--k', type=int, default=40)
    args = parser.parse_args()
    k = args.k

    ref = count_kmers_with_position(
        SeqIO.to_dict(SeqIO.parse(args.ref, "fasta")), k)
    asm = count_kmers_with_position(
        SeqIO.to_dict(SeqIO.parse(args.asm, "fasta")), k)

    n_ref = len(ref.keys())
    n_asm = len(asm.keys())
    print('n_ref', n_ref)
    print('n_asm', n_asm)

    n_both = len(set(ref.keys()) | set(asm.keys()))
    print('n_both', n_both)

    n_unique_kmers, n_unique_kmers_contained = unique_kmers_in_ref(ref, asm)
    print('n_unique', n_unique_kmers, n_unique_kmers_contained,
          n_unique_kmers_contained / n_unique_kmers)

    n_ref_only, n_asm_only = disjoint_kmers(ref, asm)
    print('n_ref_only', n_ref_only, n_ref_only / n_ref)
    print('n_asm_only', n_asm_only, n_asm_only / n_asm)

    n_same, n_different = same_count(ref, asm)
    print('n_same', n_same, n_same / n_both)
    print('n_different', n_different, n_different / n_both)


if __name__ == '__main__':
    main()
