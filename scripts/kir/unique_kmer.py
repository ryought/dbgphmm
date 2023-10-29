#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
% seqkit stat kir/hifi/10x/euler.fa
file                   format  type  num_seqs  sum_len  min_len  avg_len  max_len
kir/hifi/10x/euler.fa  FASTA   DNA          2  462,536  207,239  231,268  255,297
% seqkit stat kir/hifi/hifiasm/10x.fa
file                     format  type  num_seqs  sum_len  min_len    avg_len  max_len
kir/hifi/hifiasm/10x.fa  FASTA   DNA          3  452,671   38,370  150,890.3  207,248

% python scripts/kir/unique_kmer.py kir/hifi/t2t/KIR.fa kir/hifi/10x/euler.fa --k 40
21198 20084 0.9474478724407963
% python scripts/kir/unique_kmer.py kir/hifi/t2t/KIR.fa kir/hifi/10x/v0.final.fa --k 40
21198 19986 0.9428247947919615
% python scripts/kir/unique_kmer.py kir/hifi/t2t/KIR.fa kir/hifi/hifiasm/10x.fa --k 40
21198 15960 0.7529012170959525
% python scripts/kir/unique_kmer.py kir/hifi/t2t/KIR.fa kir/hifi/LJA/10x/10x.fa --k 40
21198 13042 0.6152467213888103
"""

from Bio import SeqIO
import argparse
from pathlib import Path
import sys
from collections import defaultdict


def disjoint_kmers(ref, asm):
    n_ref_only = 0
    n_asm_only = 0
    for kmer, count in ref.items():
        if kmer not in asm:
            n_ref_only += 1
    for kmer, count in asm.items():
        if kmer not in ref:
            n_asm_only += 1
    return n_ref_only, n_asm_only


def unique_kmers_in_ref(ref, asm):
    n_unique_kmers = 0
    n_unique_kmers_contained = 0
    for kmer, count in ref.items():
        if count == 1:
            n_unique_kmers += 1
            if kmer in asm:
                n_unique_kmers_contained += 1
    return n_unique_kmers, n_unique_kmers_contained


def count_kmers(seqs, k):
    kmers = defaultdict(int)
    for seq in seqs.values():
        s = sanitize(str(seq.seq))
        for i in range(len(s) - k + 1):
            kmer = s[i:i+k]
            kmers[kmer] += 1
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

    ref = count_kmers(SeqIO.to_dict(SeqIO.parse(args.ref, "fasta")), k)
    asm = count_kmers(SeqIO.to_dict(SeqIO.parse(args.asm, "fasta")), k)

    n_ref = len(ref.keys())
    n_asm = len(asm.keys())
    print('n_ref', n_ref)
    print('n_asm', n_asm)

    n_unique_kmers, n_unique_kmers_contained = unique_kmers_in_ref(ref, asm)
    print(n_unique_kmers, n_unique_kmers_contained,
          n_unique_kmers_contained / n_unique_kmers)

    n_ref_only, n_asm_only = disjoint_kmers(ref, asm)
    print('n_ref_only', n_ref_only, n_ref_only / n_ref)
    print('n_asm_only', n_asm_only, n_asm_only / n_asm)


if __name__ == '__main__':
    main()
