#!/bin/bash
set -Ceuo pipefail
cargo build --release
ARG="-c 20 -l 100 -p 0.01 --k-init 40 --k-final 1100 -U 10000 -N 10 -E 50 -P 2 -D 0.0 -H 0.001 -s 1 --p-infer 0.001 --use-homo-ends --start-from-true-dbg --use-true-end-nodes -m 50 -d 3"
OUTPUT="LargeK"
pjsub -x "ARG=$ARG,OUTPUT=$OUTPUT" -N $OUTPUT -o $OUTPUT.log -j sample_posterior.sh
