#!/bin/bash

cd data/em-freq

# dbgphmm -k 8 generate -l 50 > g1.fa
# dbgphmm -k 8 stat g1.fa > g1.json
# dbgphmm -k 30 sample g1.fa --length 1000 --n-reads 10 --start-from-head > g1_r1.fa
# dbgphmm -k 8 benchmark g1_r1.fa g1.fa -V 10 --init-state read-count full-em --depth-scheduler linear-gradient > g1_k8.tsv
# dbgphmm -k 16 benchmark g1_r1.fa g1.fa -V 10 --init-state read-count full-em --depth-scheduler linear-gradient > g1_k16.tsv
# dbgphmm -k 32 benchmark g1_r1.fa g1.fa -V 10 --init-state read-count full-em --depth-scheduler linear-gradient > g1_k32.tsv

# dbgphmm -k 8 generate -l 200 > g2.fa
# dbgphmm -k 8 stat g2.fa > g2.json
# dbgphmm -k 120 sample g2.fa --length 1000 --n-reads 10 --start-from-head > g2_r1.fa
# dbgphmm -k 8 benchmark g2_r1.fa g2.fa -V 10 --init-state read-count full-em --depth-scheduler linear-gradient > g2_k8.tsv
# dbgphmm -k 16 benchmark g2_r1.fa g2.fa -V 10 --init-state read-count full-em --depth-scheduler linear-gradient > g2_k16.tsv
# dbgphmm -k 32 benchmark g2_r1.fa g2.fa -V 10 --init-state read-count full-em --depth-scheduler linear-gradient > g2_k32.tsv

# dbgphmm -k 8 generate -l 500 > g3.fa
# dbgphmm -k 8 stat g3.fa > g3.json
# dbgphmm -k 250 sample g3.fa --length 1000 --n-reads 10 --start-from-head > g3_r1.fa

dbgphmm -k 16 benchmark g3_r1.fa g3.fa -V 10 --init-state read-count full-em --depth-scheduler linear-gradient > g3_k16.tsv
