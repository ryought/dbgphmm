#!/bin/sh
#PJM -g gg57
#PJM -L node=1
#PJM -L rscgrp=gg57
#PJM -L elapse=120:00:00

# or

# to submit
# pjsub -N KIR -j scripts/kir/run.sh


# if cluter
cd /work/00/gg57/j29006/dbgphmm
module load python/3.7.3

export OMP_NUM_THREADS=1

# (1) kir/haps
# 20x
# ./target/release/draft -k 40 --genome-fasta kir/haps01.fa -C 20 -L 20000 -p 0.001 -U 10 -N 10 -E 10 -H 0.0 --H0 0.0 -P 2 --output-prefix kir/haps/v0
# ./target/release/infer -k 40 -K 20000 -p 0.00001 -e 0.001 -I 50 -s 10000 --dataset-json kir/haps/v0.json --output-prefix kir/haps/v0
# 10x draft was wrong
# ./target/release/draft -k 40 --genome-fasta kir/haps/haps01.fa -C 10 -L 20000 -p 0.001 -U 10 -N 10 -E 10 -H 0.0 --H0 0.0 -P 2 --output-prefix kir/haps/10x/v0
# ./target/release/infer -k 40 -K 20000 -p 0.00001 -e 0.001 -I 50 -s 10000 --dataset-json kir/haps/10x/v0.json --output-prefix kir/haps/10x/v0
# 10x
# ./target/release/draft -k 32 -M 3 --genome-fasta kir/haps01.fa -C 10 -L 20000 -p 0.001 -U 10 -N 10 -E 10 -H 0.0 --H0 0.0 -P 2 --output-prefix kir/haps/10x/v0
# ./target/release/infer --dbg kir/haps/10x/v0.dbg -K 20000 -p 0.00001 -e 0.001 -I 50 -s 10000 --dataset-json kir/haps/10x/v0.json --output-prefix kir/haps/10x/v0

./target/release/draft -k 40 --genome-fasta kir/haps01.fa -C 10 -L 20000 --read-seed 5 -p 0.001 -m 2 -M 4 -U 0 -N 0 -E 0 -H 0.0 --H0 0.0 -P 2 --output-prefix kir/haps/10x/v0
./target/release/infer --dbg kir/haps/10x/v0.dbg -K 20000 -p 0.00001 -e 0.001 -I 50 -s 10000 --dataset-json kir/haps/10x/v0.json --output-prefix kir/haps/10x/v0

# (2) kir/hifi
# ./target/release/dbgphmm draft -k 40 -G 360000 -m 2 -M 10 -p 0.001 -P 2 -d kir/hifi/trimed.dbg -g kir/hifi/trimed.gfa kir/hifi/trimed.fa
# full
# ./target/release/dbgphmm infer -K 20000 -G 360000 -S 10000 -p 0.001 -e 0.00001 -d kir/hifi/trimed.dbg -o kir/hifi/v0 kir/hifi/trimed.fa
# 10x
# ./target/release/dbgphmm infer -K 20000 -G 360000 -S 10000 -p 0.001 -e 0.00001 -d kir/hifi/10x/v0.dbg -o kir/hifi/10x/v0 kir/hifi/10x.fa
# 10x continued
# ./target/release/dbgphmm infer -K 20000 -G 360000 -S 10000 -p 0.001 -e 0.00001 -d kir/hifi/10x/v0.k173.dbg -o kir/hifi/10x/v0 kir/hifi/10x.fa
# 20x
# ./target/release/dbgphmm infer -K 20000 -G 360000 -S 10000 -p 0.001 -e 0.00001 -d kir/hifi/20x/v0.dbg -o kir/hifi/20x/v0 kir/hifi/20x.fa
# 20x continued
# ./target/release/dbgphmm infer -K 20000 -G 360000 -S 10000 -p 0.001 -e 0.00001 -d kir/hifi/20x/v0.k53.dbg -o kir/hifi/20x/v0 kir/hifi/20x.fa
./target/release/dbgphmm infer -K 20000 -G 360000 -S 10000 -p 0.001 -e 0.00001 -d kir/hifi/20x/v0.k81.dbg -o kir/hifi/20x/v0 kir/hifi/20x.fa
