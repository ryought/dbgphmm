#!/bin/bash
set -Ceuo pipefail
cargo build --release

for U in 10 20 50 100 200 500
do
  N=$((1000 / $U))
  # echo "U=$U N=$N"
  for H in 0.0 0.001 0.01 0.05 0.1
  do
    for S in 2 3
    do
      ARG="-c 20 -l 100 -p 0.01 --k-init 12 --k-final 15 -U $U -N $N -E 50 -P 2 -D 0.0 -H $H --sigma 100 -d 10 -m 50 --use-true-end-nodes --start-from-true --use-true-dbg -s $S"
      OUTPUT="U${U}H${H//./}S${S}"
      echo $ARG $OUTPUT
      pjsub -x "ARG=$ARG,OUTPUT=$OUTPUT" -N $OUTPUT -o $OUTPUT.log -j sample_posterior.sh
    done
  done
done
