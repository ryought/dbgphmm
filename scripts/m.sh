#!/bin/sh
#PJM -g gg57
#PJM -L node=1
#PJM -L elapse=24:00:00
#PJM -L rscgrp=gg57

if [ $1 = "sub" ]
then
  # submitter
  echo "pjsub"

  cargo build --release
  pjsub -N p01_u1k_n10 -j scripts/m.sh
else
  #
  # runner
  #
  cd /work/00/gg57/j29006/dbgphmm
  module load python/3.7.3
  ./target/release/draft -k 40 -C 20 -L 2000 -p 0.001 -U 1000 -N 10 -E 300 -H 0.01 -P 2 --output-prefix m/p01_u1k_n10
  ./target/release/infer -k 40 -K 2000 -p 0.0001 -e 0.001 -I 50 --dataset-json m/p01_u1k_n10.json --output-prefix m/p01_u1k_n10
fi
