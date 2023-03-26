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
  for N in 1 2 3 4 5 6
  do
    echo $N
    # pjsub -x N=$N -N N$N -j -o n/N$N scripts/n.sh
    pjsub -x N=$N -N piN$N -j -o n/piN$N scripts/n.sh
  done
else
  #
  # runner
  #
  cd /work/00/gg57/j29006/dbgphmm
  module load python/3.7.3
  echo "Running N=$N"
  # ./target/release/draft -k 20 -C 20 -L 1000 -p 0.001 -U 500 -N $N -E 300 -H 0.02 -P 2 --output-prefix n/p01_u500_n$N
  # ./target/release/infer -k 20 -K 1000 -p 0.0001 -I 50 --dataset-json n/p01_u500_n$N.json --output-prefix n/p01_u500_n$N
  ./target/release/infer -k 20 -K 1000 -p 0.00001 -I 50 --dataset-json n/p01_u500_n$N.json --output-prefix n/p01_u500_n${N}_pi
fi
