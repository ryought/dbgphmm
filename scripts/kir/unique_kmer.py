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


def has_kmer(asm, kmer):
    has = False
    for hap in asm.values():
        if kmer in hap.seq:
            has = True

    return has


def main():
    parser = argparse.ArgumentParser(description='unique k-mer analysis')
    parser.add_argument('ref', type=Path, help='FASTA')
    parser.add_argument('asm', type=Path, help='FASTA')
    parser.add_argument('--k', type=int, default=100)
    args = parser.parse_args()
    k = args.k

    ref = SeqIO.to_dict(SeqIO.parse(args.ref, "fasta"))
    asm = SeqIO.to_dict(SeqIO.parse(args.asm, "fasta"))

    kmers = defaultdict(int)
    for hap in ref.values():
        for i in range(len(hap) - k + 1):
            kmer = str(hap[i:i+k].seq)
            kmers[kmer] += 1

    n_unique_kmers = 0
    n_unique_kmers_contained = 0
    for kmer, count in kmers.items():
        if count == 1:
            n_unique_kmers += 1
            has = has_kmer(asm, kmer)
            if has:
                n_unique_kmers_contained += 1

    print(n_unique_kmers, n_unique_kmers_contained,
          n_unique_kmers_contained / n_unique_kmers)


if __name__ == '__main__':
    main()
