#!/bin/bash
set -Ceuo pipefail
cargo build --release

U=20
N=50
H=0.001
for S in 0 1 2
do
  p_infer=0.001
  # raw
  # ARG="-c 20 -l 100 -p 0.01 --k-init 12 --k-final 100 -U $U -N $N -E 50 -P 2 -D 0.0 -H $H --sigma 100 -d 10 -m 50 --use-true-end-nodes -s $S --p-infer $p_infer --p0 0.8 --use-homo-ends"
  # OUTPUT="LowPInfer_U${U}N${N}H${H//./}S${S}PInfer${p_infer}"
  # pjsub -x "ARG=$ARG,OUTPUT=$OUTPUT" -N $OUTPUT -o $OUTPUT.log -j sample_posterior.sh

  # fromtrue
  # ARG="-c 20 -l 100 -p 0.01 --k-init 12 --k-final 100 -U $U -N $N -E 50 -P 2 -D 0.0 -H $H --sigma 100 -d 10 -m 50 --use-true-end-nodes -s $S --p-infer $p_infer --p0 0.8 --use-homo-ends --start-from-true"
  # OUTPUT="LowPInfer_U${U}N${N}H${H//./}S${S}FromTruePInfer${p_infer}"
  # pjsub -x "ARG=$ARG,OUTPUT=$OUTPUT" -N $OUTPUT -o $OUTPUT.log -j sample_posterior.sh
done

for U in 1000 500 100 50 20
do
  L=1000
  N=$(($L / $U))
  for S in 0 1 2
  do
    for p_infer in 0.005 0.001
    do
      ARG="-c 20 -l 100 -p 0.01 --k-init 12 --k-final 900 -U $U -N $N -E 50 -P 2 -D 0.0 -H $H --sigma 100 -d 10 -m 50 --use-true-end-nodes -s $S --p-infer $p_infer --p0 0.8 --use-homo-ends --start-from-true-dbg"
      OUTPUT="LowPInfer_U${U}N${N}H${H//./}S${S}FromTrueDbgPInfer${p_infer}"
      pjsub -x "ARG=$ARG,OUTPUT=$OUTPUT" -N $OUTPUT -o $OUTPUT.log -j sample_posterior.sh
    done
  done
done

# L=1000
# for U in 10 20 50 100 200 500
# do
#   N=$(($L / $U))
#   # echo "U=$U N=$N"
#   for H in 0.0 0.001 0.01 0.05 0.1
#   do
#     for S in 0 1 2
#     do
#       ARG="-c 20 -l 100 -p 0.01 --k-init 12 --k-final 100 -U $U -N $N -E 50 -P 2 -D 0.0 -H $H --sigma 100 -d 10 -m 50 --use-true-end-nodes -s $S"
#       OUTPUT="LowPInfer_U${U}N${N}H${H//./}S${S}"
#       echo $ARG $OUTPUT
#       pjsub -x "ARG=$ARG,OUTPUT=$OUTPUT" -N $OUTPUT -o $OUTPUT.log -j sample_posterior.sh
#     done
#   done
# done
