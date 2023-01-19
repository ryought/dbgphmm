#!/bin/bash
set -Ceuo pipefail
cargo build --release
# cargo run --release --bin batch -- -H 0.01 -U 100 -N 10 -P 2
for U in 1000 500 100 50 20
do
  L=1000
  N=$(($L / $U))
  for H in 0.001 0.005
  do
    for p in 0.001 0.005 0.01
    do
      ARG="-U $U -N $N -H $H -P 2 -p $p --try-all-k"
      OUTPUT="U${U}N${N}H${H//./}P2p${p//./}"
      echo $ARG $OUTPUT
      pjsub -x "ARG=$ARG,OUTPUT=$OUTPUT" -N $OUTPUT -o $OUTPUT.log -j job.sh
    done
  done
done
