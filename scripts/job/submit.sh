#!/bin/bash
set -Ceuo pipefail
case $1 in
  "case1")
    # good, easiest case
    ARG="-c 10 -l 100 --k-init 32 --k-final 100 -U 1000 -N 1 -E 50 -P 2 -H 0.01 --sigma 100 -m 10 --use-fragment-read";;
  "case2")
    # good
    ARG="-c 20 -l 100 -p 0.01 --k-init 16 --k-final 100 -U 200 -N 5 -E 50 -P 1 -D 0.01 -H 0.01 --sigma 100 -m 10 --use-fragment-read";;
  "case2t")
    # good
    ARG="-c 20 -l 100 -p 0.01 --k-init 16 --k-final 100 -U 200 -N 5 -E 50 -P 1 -D 0.01 -H 0.01 --sigma 100 -m 10 --use-fragment-read --start-from-true";;
  "case3k12")
    # fail in k=12
    # use in sampling test stopped local minimum
    # ARG="-c 20 -l 100 -p 0.01 --k-init 12 --k-final 14 -U 50 -N 20 -E 50 -P 2 -D 0.05 -H 0.05 --sigma 100 -d 5 -m 40 --use-fragment-read";;
    ARG="-c 20 -l 100 -p 0.01 --k-init 12 --k-final 14 -U 50 -N 20 -E 50 -P 2 -D 0.05 -H 0.05 --sigma 100 -d 10 -m 50 --use-fragment-read --dbgviz-output case3k12.json";;
  "case3k12t")
    # fail in k=13
    ARG="-c 20 -l 100 -p 0.01 --k-init 12 --k-final 100 -U 50 -N 20 -E 50 -P 2 -D 0.05 -H 0.05 --sigma 100 -m 10 --use-fragment-read --start-from-true";;
  # k=16
  # missing kmer
  "case3k16")
    ARG="-c 20 -l 100 -p 0.01 --k-init 16 --k-final 100 -U 50 -N 20 -E 50 -P 2 -D 0.05 -H 0.05 --sigma 100 -m 10 --use-fragment-read";;
  "case3k16t")
    ARG="-c 20 -l 100 -p 0.01 --k-init 16 --k-final 100 -U 50 -N 20 -E 50 -P 2 -D 0.05 -H 0.05 --sigma 100 -m 10 --use-fragment-read --start-from-true";;
  "case3L")
    # deprecated same as case3k16t
    ARG="-c 20 -l 100 -p 0.01 --k-init 16 --k-final 100 -U 50 -N 20 -E 50 -P 1 -D 0.05 -H 0.05 --sigma 100 -m 10 --use-fragment-read";;
  "case4")
    ARG="-c 15 -l 100 -p 0.001 --k-init 16 --k-final 100 -U 50 -N 10 -E 50 -P 2 -D 0.05 -H 0.05 --sigma 100 -m 10 --use-fragment-read";;
  "case4t")
    ARG="-c 15 -l 100 -p 0.001 --k-init 16 --k-final 100 -U 50 -N 10 -E 50 -P 2 -D 0.05 -H 0.05 --sigma 100 -m 10 --use-fragment-read --start-from-true";;
  # "case4a")
  #   # missing kmer
  #   ARG="-c 15 -l 100 -p 0.001 --k-init 24 --k-final 100 -U 50 -N 10 -E 50 -P 2 -D 0.05 -H 0.05 --sigma 100 -m 3 --use-fragment-read";;
  # "case4b")
  #   # missing kmer
  #   ARG="-c 15 -l 100 -p 0.001 --k-init 32 --k-final 100 -U 50 -N 10 -E 50 -P 2 -D 0.05 -H 0.05 --sigma 100 -m 3 --use-fragment-read";;
  # "case4c")
  #   # missing
  #   ARG="-c 15 -l 100 -p 0.001 --k-init 12 --k-final 100 -U 50 -N 10 -E 50 -P 2 -D 0.05 -H 0.05 --sigma 100 -m 3 --use-fragment-read";;
  "case5")
    # good easy case
    ARG="-c 10 -l 200 -p 0.001 --k-init 40 --k-final 200 -U 2000 -N 1 -E 0 -P 2 -D 0.05 -H 0.05 --sigma 100 -m 3 --use-fragment-read";;
  "case6")
    # 2000bp diploid + 1% 200bp 20x read
    # k=16 some k-mer is missing
    # k=12 still running.. too many neighbors
    ARG="-c 20 -l 200 -p 0.01 --k-init 12 --k-final 220 -U 2000 -N 1 -E 0 -P 2 -D 0.05 -H 0.05 --sigma 100 -m 10 --use-fragment-read";;
  "case6F")
    # 2000bp diploid + 1% 2000bp 20x FullLength read
    # fail k=16 missing kmer
    ARG="-c 20 -l 2000 -p 0.01 --k-init 16 --k-final 2000 -U 2000 -N 1 -E 0 -P 2 -D 0.05 -H 0.05 --sigma 100 -m 10";;
  "case7")
    # 500bp x4 diploid + 1% 200bp 20x read
    # still running
    ARG="-c 20 -l 200 -p 0.01 --k-init 16 --k-final 220 -U 500 -N 4 -E 50 -P 2 -D 0.05 -H 0.05 --sigma 100 -m 10 --use-fragment-read";;
  "case7F")
    # 500bp x4 diploid + 1% 200bp 20x FullLength read
    # fail in k=16 purge, incomplete search
    ARG="-c 20 -l 2000 -p 0.01 --k-init 16 --k-final 2000 -U 500 -N 4 -E 50 -P 2 -D 0.05 -H 0.05 --sigma 100 -m 10";;
  "case8")
    # 100bp x20 haploid + 1% 200bp 20x read
    # k=16 some kmer is missing
    # k=12 is also missing!
    # k=8 is missing in first run
    ARG="-c 20 -l 200 -p 0.01 --k-init 8 --k-final 200 -U 100 -N 20 -E 50 -P 1 -D 0.01 -H 0.01 --sigma 100 -m 10 --use-fragment-read";;
  "case9")
    # accurate read case
    # 100bp x20 diploid + 0.1% 500bp 10x read
    # good? but takes long time
    ARG="-c 10 -l 500 -p 0.001 --k-init 32 --k-final 500 -U 100 -N 20 -E 50 -P 2 -D 0.01 -H 0.01 --sigma 100 -m 10 --use-fragment-read";;
  "case9t")
    ARG="-c 10 -l 500 -p 0.001 --k-init 32 --k-final 500 -U 100 -N 20 -E 50 -P 2 -D 0.01 -H 0.01 --sigma 100 -m 10 --use-fragment-read --start-from-true";;
  *)
    echo "unknown case id" && exit 1;;
esac
OUTPUT=$1
echo $ARG
echo $OUTPUT
cargo build --release
pjsub -x "ARG=$ARG,OUTPUT=$OUTPUT" -N $OUTPUT -o $OUTPUT.log -j sample_posterior.sh
# pjsub -N $OUTPUT -o $OUTPUT.log -j sample_posterior.sh
