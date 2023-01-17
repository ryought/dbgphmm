#!/bin/bash
set -Ceuo pipefail

for U in 10 20 50 100 200 500
do
  L=500
  N=$(($L / $U))
  for H in 0.0 0.001 0.01 0.05 0.1
  do
    for S in 0 1 2 3 4
    do
      OUTPUT="U${U}N${N}H${H//./}S${S}"

      # .kmers
      # echo $OUTPUT
      # awk "\$1==\"K\" && \$5==1 && \$9 > 0.5" $OUTPUT

      # .summary
      for K in $(seq 12 15)
      do
        echo $OUTPUT $K $(awk "\$1==\"N\" && \$2==$K" $OUTPUT | tail -n 1)
        # echo $OUTPUT $K $(awk "\$1==\"K\" && \$2==$K" $OUTPUT | wc -l)
      done
    done
  done
done
