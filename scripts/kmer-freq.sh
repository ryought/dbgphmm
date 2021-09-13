#!/bin/bash
# count-base
# good?

cd data/kmer-freq

dbgphmm -k 8 generate -l 50 > g1.fa
dbgphmm -k 8 stat g1.fa > g1.json
dbgphmm -k 30 sample g1.fa --length 1000 --n-reads 10 --start-from-head > g1_r1.fa
dbgphmm -k 8 benchmark g1_r1.fa g1.fa -V 10 --init-state zero freq-em > g1_k8.tsv
dbgphmm -k 16 benchmark g1_r1.fa g1.fa -V 10 --init-state zero freq-em > g1_k16.tsv
dbgphmm -k 32 benchmark g1_r1.fa g1.fa -V 10 --init-state zero freq-em > g1_k32.tsv

dbgphmm -k 8 generate -l 200 > g2.fa
dbgphmm -k 8 stat g2.fa > g2.json
dbgphmm -k 120 sample g2.fa --length 1000 --n-reads 10 --start-from-head > g2_r1.fa
dbgphmm -k 8 benchmark g2_r1.fa g2.fa -V 10 --init-state zero freq-em > g2_k8.tsv
dbgphmm -k 16 benchmark g2_r1.fa g2.fa -V 10 --init-state zero freq-em > g2_k16.tsv
dbgphmm -k 32 benchmark g2_r1.fa g2.fa -V 10 --init-state zero freq-em > g2_k32.tsv
