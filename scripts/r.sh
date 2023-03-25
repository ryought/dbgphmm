#!/bin/sh
#PJM -g gg57
#PJM -L node=1
#PJM -L elapse=24:00:00
#PJM -L rscgrp=gg57

#
# To submit:
# pjsub -X KEY=VALUE,KEY=VALUE -N ${JOB_NAME} -o ${LOG_FILE} -j --comment hoge job.sh
#
cd /work/00/gg57/j29006/dbgphmm
module load python/3.7.3

# ./target/release/draft -k 16 -C 20 -L 500 -p 0.01 -U 1000 -N 1 -E 200 -H 0.03 -P 2 --output-prefix r/u1k_n1
# ./target/release/infer -k 16 -K 500 -p 0.001 --dataset-json r/u1k_n1.json --output-prefix r/u1k_n1 -I 50

# ./target/release/draft -k 20 -C 20 -L 400 -p 0.005 -U 200 -N 10 -E 200 -H 0.03 -P 2 --output-prefix r/u200_n10
# ./target/release/infer -k 20 -K 40 -p 0.001 --dataset-json r/u200_n10.json --output-prefix r/u200_n10_v3 -I 50
# ./target/release/infer -k 20 -K 40 -p 0.0001 --dataset-json r/u200_n10.json --output-prefix r/u200_n10_p00001 -I 50

KEY="u10k_n2"
./target/release/draft -k 20 -C 20 -L 20000 -p 0.005 -U 10000 -N 2 -E 1000 -H 0.03 -P 2 --output-prefix r/$KEY
./target/release/infer -K 20000 -p 0.001 -I 50 --dataset-json r/$KEY.json --dbg r/$KEY.dbg --output-prefix r/$KEY
