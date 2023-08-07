#!/bin/sh
#PJM -g gg57
#PJM -L rscgrp=short
#PJM -L node=1
#PJM -L elapse=8:00:00

cd /work/00/gg57/j29006/dbgphmm
module load python/3.7.3
export OMP_NUM_THREADS=1

# ./target/release/sample -p 0.00001 -s 10000 --dbg L_inspect/v0.k44.dbg --dataset-json L_inspect/v0.json --edges 571,2770,2889
./target/release/sample -p 0.00001 -s 10000 --dbg L_inspect/v0.k44.dbg --dataset-json L/v0.json --edges 571 --inspect-filename L_inspect/v0.k44.e571.allowk300.inspect --k 300 --allow-zero-edge
./target/release/sample -p 0.00001 -s 10000 --dbg L_inspect/v0.k44.dbg --dataset-json L/v0.json --edges 571 --inspect-filename L_inspect/v0.k44.e571.disallowk300.inspect --k 300

# ./target/release/mapping --dbg L/v0.k44.dbg --dataset L/v0.json -e 0.00001 -a 40 -r 30 --phmm normal --map-input L/v0.k44.map > L/v0.k44.bmap
# ./target/release/mapping --dbg L/v0.k44.dbg --dataset L/v0.json -e 0.00001 -a 40 -r 30 --phmm normal --map-input L/v0.k44.map > L/v0.k44.bmap
