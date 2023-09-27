#!/bin/sh
#PJM -g gg57
#PJM -L rscgrp=regular
#PJM -L node=1
#PJM -L elapse=48:00:00

cd /work/00/gg57/j29006/dbgphmm
module load python/3.7.3
export OMP_NUM_THREADS=1
# ./target/release/draft -k 40 -C 20 -L 10000 -p 0.001 -U 2000 -N 50 -E 10000 -H 0.001 --H0 0.001 -P 2 --output-prefix t/t2
./target/release/draft -k 40 -C 20 -L 10000 -p 0.001 -U 2000 -N 20 -E 5000  -H 0.001 --H0 0.001 -P 2 --output-prefix t/t8
./target/release/infer -k 40 -p 0.00001 -K 10000 -e 0.001 -I 50 -s 10000 --dbg t/t8.dbg --dataset-json t/t8.json --output-prefix t/t8
