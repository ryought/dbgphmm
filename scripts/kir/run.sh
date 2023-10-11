#!/bin/sh
#PJM -g gg57
#PJM -L node=1
#PJM -L rscgrp=gg57
#PJM -L elapse=120:00:00

# to submit
# pjsub -N KIR -j scripts/kir/run.sh

# if cluter
cd /work/00/gg57/j29006/dbgphmm
module load python/3.7.3

export OMP_NUM_THREADS=1

# (1) kir/haps 10x
./target/release/draft -k 40 --genome-fasta kir/haps01.fa -C 10 -L 20000 --read-seed 5 -p 0.001 -m 2 -M 4 -U 0 -N 0 -E 0 -H 0.0 --H0 0.0 -P 2 --output-prefix kir/haps/10x/v0
./target/release/infer --dbg kir/haps/10x/v0.dbg -K 20000 -p 0.00001 -e 0.001 -I 50 -s 10000 --dataset-json kir/haps/10x/v0.json --output-prefix kir/haps/10x/v0

# (2) kir/hifi 10x
./target/release/dbgphmm draft -k 40 -G 360000 -m 2 -M 5 -p 0.001 -P 2 -d kir/hifi/10x/v0.dbg -g kir/hifi/10x/v0.gfa kir/hifi/10x.fa
./target/release/dbgphmm infer -K 20000 -G 360000 -S 10000 -p 0.0003 -e 0.0003 --p0 0.99 -d kir/hifi/10x/v0.dbg -o kir/hifi/10x/v0 kir/hifi/10x.fa
# 20x
# ./target/release/dbgphmm infer -K 20000 -G 360000 -S 10000 -p 0.001 -e 0.00001 -d kir/hifi/20x/v0.dbg -o kir/hifi/20x/v0 kir/hifi/20x.fa
# ./target/release/dbgphmm draft -k 40 -G 360000 -m 2 -M 8 -p 0.001 -P 2 -d kir/hifi/20x/v0.dbg -g kir/hifi/20x/v0.gfa kir/hifi/20x.fa
