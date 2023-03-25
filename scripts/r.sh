#!/bin/sh
#PJM -g gg57
#PJM -L node=1
#PJM -L elapse=24:00:00
#PJM -L rscgrp=gg57

#
# To submit:
#
#
cd /work/00/gg57/j29006/dbgphmm
module load python/3.7.3
./target/release/draft -k 16 -C 20 -L 500 -p 0.01 -U 1000 -N 1 -E 200 -H 0.03 -P 2 --output-prefix r/u1k_n1
sleep 5
./target/release/infer -k 16 -K 500 -p 0.001 --dataset-json r/u1k_n1.json --output-prefix r/u1k_n1 -I 50

# ./target/release/draft -k 20 -C 20 -L 400 -p 0.005 -U 200 -N 10 -E 200 -H 0.03 -P 2 --output-prefix r/repeat_u200
# ./target/release/infer -k 20 -K 40 -p 0.001 --dataset-json r/repeat_u200.json --output-prefix r/repeat_u200_v2 -I 50
