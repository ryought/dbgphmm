#!/bin/sh
#PJM -g gg57
#PJM -L rscgrp=gg57
#PJM -L node=1
#PJM -L elapse=48:00:00

if [ $1 = "sub" ]
then
  # submitter
  echo "pjsub"

  cargo build --release --features intel --no-default-features
  mkdir -p L
  # pjsub -x H=$H -N L -j scripts/l.sh
  pjsub -N L -j scripts/l.sh
else
  #
  # runner
  #
  cd /work/00/gg57/j29006/dbgphmm
  module load python/3.7.3
  export OMP_NUM_THREADS=1
  echo "Running L"
  # 10kb unit repeated 10 times
  # 20x 20kb read
  K=20000
  ./target/release/draft -k 40 -C 20 -L $K -p 0.001 -U 10000 -N 10 -E 500 -H 0.02 --H0 0.02 -P 2 --output-prefix L/v0
  ./target/release/infer -k 40 -K 100 -p 0.00001 -e 0.001 -I 50 -s 10000 --km $K --kr $K --dataset-json L/v0.json --output-prefix L/v0
fi
