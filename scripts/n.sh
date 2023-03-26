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
  mkdir -p n
  pjsub -x N=1 -N N1 -j scripts/n.sh
else
  # runner
  #
  cd /work/00/gg57/j29006/dbgphmm
  module load python/3.7.3
  echo "Running N=$N"
  ./target/release/draft -k 20 -C 20 -L 1000 -p 0.01 -U 500 -N $N -E 300 -H 0.02 -P 2 --output-prefix n/u500_n$N
  ./target/release/infer -k 20 -K 1000 -p 0.001 -I 50 --dataset-json n/u500_n$N.json --output-prefix n/u500_n$N
fi

