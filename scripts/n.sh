#!/bin/sh
#PJM -g gg57
#PJM -L rscgrp=short
#PJM -L node=1
#PJM -L elapse=8:00:00

#
# rscgrp=gg57 or short

if [ $1 = "sub" ]
then
  # submitter
  echo "pjsub"

  cargo build --release --features intel --no-default-features
  mkdir -p n3

  if [ $2 = "u100" ]
  then
    # u100
    for N in 10 20 50 100
    do
      echo $N
      pjsub -x N=$N,U=100 -N n3_p01_u100_n$N -j scripts/n.sh
    done
  elif [ $2 = "u500" ]
  then
    # u500
    for N in 1 2 3 4 5 6 7 8 9 10 11 12
    do
      echo $N
      pjsub -x N=$N,U=500 -N n3_p01_u500_n$N -j scripts/n.sh
    done
  elif [ $2 = "u20" ]
  then
    # u20
    for N in 50 100 200
    do
      echo $N
      pjsub -x N=$N,U=20 -N n3_p01_u20_n$N -j scripts/n.sh
    done
  fi
else
  #
  # runner
  #
  cd /work/00/gg57/j29006/dbgphmm
  module load python/3.7.3
  export OMP_NUM_THREADS=1
  echo "Running N=$N U=$U"

  # p01
  if [ $U = 500 ]
  then
    # u500
    ./target/release/draft -k 40 -C 20 -L 1000 -p 0.001 -U 500 -N $N -E 300 -H 0.02 --H0 0.02 -P 2 --output-prefix n3/p01_u500_n$N
    ./target/release/infer -k 40 -K 1000 -p 0.0001 -e 0.001 -I 50 -s 500 --dataset-json n3/p01_u500_n$N.json --output-prefix n3/p01_u500_n$N
  elif [ $U = 20 ]
  then
    # u20
    ./target/release/draft -k 40 -C 20 -L 1000 -p 0.001 -U 20 -N $N -E 300 -H 0.02 --H0 0.02 -P 2 --output-prefix n3/p01_u20_n$N
    ./target/release/infer -k 40 -K 1000 -p 0.00001 -e 0.001 -I 50 -s 500 --dataset-json n3/p01_u20_n$N.json --output-prefix n3/p01_u20_n$N
  elif [ $U = 100 ]
  then
    # u100
    ./target/release/draft -k 40 -C 20 -L 1000 -p 0.001 -U 100 -N $N -E 300 -H 0.02 --H0 0.02 -P 2 --output-prefix n3/p01_u100_n$N
    ./target/release/infer -k 40 -K 1000 -p 0.00001 -e 0.001 -I 50 -s 500 --dataset-json n3/p01_u100_n$N.json --output-prefix n3/p01_u100_n$N
  fi

  # p1
  # ./target/release/draft -k 16 -C 20 -L 1000 -p 0.01 -U 500 -N $N -E 300 -H 0.02 --H0 0.02 -P 2 --output-prefix n3/p1_u500_n$N
  # ./target/release/infer -k 16 -K 1000 -p 0.0001 -e 0.01 -I 100 -s 500 --km 1000 --kr 1000 --dataset-json n3/p1_u500_n$N.json --output-prefix n3/p1_u500_n$N
fi
