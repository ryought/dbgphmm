#!/bin/bash
set -Ceuo pipefail
#!
case $1 in
  "case1")
    ARG="-c 10 -l 100 --k-init 32 --k-final 100 -U 1000 -N 1 -E 50 -P 2 -H 0.01 --sigma 100 -m 3 --use-fragment-read";;
  "case2")
    ARG="-c 20 -l 100 -p 0.01 --k-init 16 --k-final 100 -U 200 -N 5 -E 50 -P 1 -D 0.01 -H 0.01 --sigma 100 -m 3 --use-fragment-read";;
  "case3")
    ARG="-c 20 -l 100 -p 0.01 --k-init 12 --k-final 100 -U 50 -N 20 -E 50 -P 2 -D 0.05 -H 0.05 --sigma 100 -m 3 --use-fragment-read";;
  "case4")
    ARG="-c 15 -l 100 -p 0.001 --k-init 16 --k-final 100 -U 50 -N 10 -E 50 -P 2 -D 0.05 -H 0.05 --sigma 100 -m 3 --use-fragment-read";;
  "case5")
    ARG="-c 10 -l 200 -p 0.001 --k-init 40 --k-final 200 -U 2000 -N 1 -E 0 -P 2 -D 0.05 -H 0.05 --sigma 100 -m 3 --use-fragment-read";;
  *)
    echo "unknown case id" && exit 1;;
esac
OUTPUT=$1
echo $ARG
echo $OUTPUT
pjsub -x "ARG=$ARG,OUTPUT=$OUTPUT" -N $OUTPUT -o $OUTPUT.log -j sample_posterior.sh
# pjsub -N $OUTPUT -o $OUTPUT.log -j sample_posterior.sh
