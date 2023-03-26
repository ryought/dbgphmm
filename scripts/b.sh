#!/bin/sh
#PJM -g gg57
#PJM -L node=1
#PJM -L elapse=24:00:00
#PJM -L rscgrp=gg57

if [ $1 = "sub" ]
then
  #
  # submitter
  #
  cargo build --release
  for N in 1 2 4
  do
    echo $N
    pjsub -x N=$N -N bn$N -j -o b/n$N scripts/b.sh
  done
else
  #
  # runner
  #
  cd /work/00/gg57/j29006/dbgphmm
  module load python/3.7.3
  echo "Running N=$N"
  # 10k
  ./target/release/draft -k 40 -C 20 -L 20000 -p 0.001 -U 10000 -N $N -E 1000 -H 0.02 -P 2 --output-prefix b/p01_u10k_n$N
  ./target/release/infer -k 40 -K 20000 -p 0.00001 -I 50 --dataset-json b/p01_u10k_n$N.json --output-prefix b/p01_u10k_n${N}
fi
