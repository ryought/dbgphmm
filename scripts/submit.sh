#!/bin/bash
set -u

for L in 50 100 500
do
  for S in 0 1 2
  do
    qsub -v L=$L -v S=$S -N "dbgphmm_L${L}_S${S}" scripts/benchmark.sh
  done
done
