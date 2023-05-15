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
  mkdir -p n
  # for N in 1 2 3 4 5 6
  # for N in 7 8 9 10 11 12
  for N in 10 11 12
  do
    echo $N
<<<<<<< Updated upstream
    # add git hash
    pjsub -x N=$N -N e_p01_n$N -j scripts/n.sh
=======
    # Name:
    # - key
    # - git-hash
    # Parameter:
    # - N
    # Compile before run
    pjsub -x N=$N -N p01_n$N -j scripts/n.sh
>>>>>>> Stashed changes
  done
else
  #
  # runner
  #
  cd /work/00/gg57/j29006/dbgphmm
  module load python/3.7.3
  export OMP_NUM_THREADS=1
  echo "Running N=$N"
  # ./target/release/draft -k 20 -C 20 -L 1000 -p 0.001 -U 500 -N $N -E 300 -H 0.02 -P 2 --output-prefix n/p01_u500_n$N
  # ./target/release/infer -k 20 -K 1000 -p 0.0001 -I 50 --dataset-json n/p01_u500_n$N.json --output-prefix n/p01_u500_n$N
  # ./target/release/infer -k 20 -K 1000 -p 0.00001 -I 50 --dataset-json n/p01_u500_n$N.json --output-prefix n/p01_u500_n${N}_pi
  # ./target/release/infer -k 20 -K 1000 -p 0.00001 -I 50 --dataset-json n/p01_u500_n$N.json --output-prefix n/p01_u500_n${N}_pi_v2
  # ./target/release/infer -k 40 -K 1000 -p 0.00001 -I 50 --dataset-json n/p01_u500_n$N.json --output-prefix n/p01_u500_n${N}_pi_v3

  # store correspondence output-data and day-or-parameters
  #

  # p01
  ./target/release/draft -k 40 -C 20 -L 1000 -p 0.001 -U 500 -N $N -E 300 -H 0.02 --H0 0.02 -P 2 --output-prefix n/p01_u500_n$N
  ./target/release/infer -k 40 -K 1000 -p 0.0001 -e 0.001 -I 50 -s 500 --dataset-json n/p01_u500_n$N.json --output-prefix n/s500_p01_u500_n$N
  # ./target/release/infer -k 40 -K 1000 -p 0.001 -e 0.001 -I 50 -s 500 --dataset-json n/p01_u500_n$N.json --output-prefix n/f_s500_p01_u500_n$N

  # p1
  # ./target/release/draft -k 16 -C 20 -L 1000 -p 0.01 -U 500 -N $N -E 300 -H 0.02 --H0 0.02 -P 2 --output-prefix n/p1_u500_n$N
  # ./target/release/infer -k 16 -K 1000 -p 0.0001 -e 0.01 -I 100 -s 500 --dataset-json n/p1_u500_n$N.json --output-prefix n/s500_p1_u500_n$N
fi
