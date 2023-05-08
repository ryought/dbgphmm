#!/bin/sh
#PJM -g gg57
#PJM -L rscgrp=regular
#PJM -L node=1
#PJM -L elapse=48:00:00

if [ $1 = "sub" ]
then
  # submitter
  echo "pjsub"
  cargo build --release --features intel --no-default-features
  mkdir -p h
  pjsub -N h -j scripts/h.sh
else
  #
  # runner
  #
  cd /work/00/gg57/j29006/dbgphmm
  module load python/3.7.3

  export OMP_NUM_THREADS=1
  ./target/release/draft -k 40 -C 20 -L 10000 -p 0.001 -U 5000 -N 10 -E 300 -H 0.02 --H0 0.02 -P 2 --output-prefix h/u5k_n10
  ./target/release/infer -k 40 -K 10000 -p 0.0001 -e 0.001 -I 50 -s 5000 --dataset-json h/u5k_n10.json --output-prefix h/u5k_n10
fi
