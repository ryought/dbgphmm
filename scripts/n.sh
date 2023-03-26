#!/bin/sh
#PJM -g gg57
#PJM -L node=1
#PJM -L elapse=24:00:00
#PJM -L rscgrp=gg57

#
# To submit:
#
# cargo build --release
# pjsub -X N=1 -N N1 -j n.sh
#

cd /work/00/gg57/j29006/dbgphmm
module load python/3.7.3
echo "Running N=$N"
./target/release/draft -k 20 -C 20 -L 1000 -p 0.01 -U 500 -N $N -E 300 -H 0.02 -P 2 --output-prefix n/u500_n$N
./target/release/infer -K 1000 -p 0.001 -I 50 --dataset-json n/u500_n$N.json --output-prefix n/u500_n$N
