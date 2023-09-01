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
  mkdir -p L1
  pjsub -N L1 -j scripts/l.sh
else
  #
  # runner
  #
  cd /work/00/gg57/j29006/dbgphmm
  module load python/3.7.3
  export OMP_NUM_THREADS=1
  echo "Running L1"
  # 10kb unit repeated 10 times
  # 20x 20kb read
  K=20000
  ./target/release/draft -k 40 -C 20 -L $K -p 0.001 -U 10000 -N 10 -E 500 -H 0.02 --H0 0.02 -P 2 --output-prefix L1/H2
  ./target/release/infer -k 40 -K $K -p 0.00001 -e 0.001 -I 50 -s 10000 --dataset-json L1/H2.json --output-prefix L1/H2
  # ./target/release/draft -k 40 -C 20 -L $K -p 0.001 -U 10000 -N 10 -E 500 -H 0.01 --H0 0.01 -P 2 --output-prefix L2/H1
  # ./target/release/infer -k 40 -K $K -p 0.00001 -e 0.001 -I 50 -s 10000 --dataset-json L2/H1.json --output-prefix L2/H1
  # ./target/release/draft -k 40 -C 20 -L $K -p 0.001 -U 10000 -N 10 -E 500 -H 0.05 --H0 0.05 -P 2 --output-prefix L3/H5
  # ./target/release/infer -k 40 -K $K -p 0.00001 -e 0.001 -I 50 -s 10000 --dataset-json L3/H5.json --output-prefix L3/H5
fi
