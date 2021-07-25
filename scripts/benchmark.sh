#!/bin/bash
set -u

DEBUG='cargo run --'
RELEASE='./target/release/dbgphmm'

function freqem () {
  LENGTH=$1
  N_HAP=$2
  SEED=$3
  ERROR_RATE=$4
  DEPTH=$5
  K=$6

  GENOME_DIR="data/L${LENGTH}_N${N_HAP}_S${SEED}"
  GENOME=$GENOME_DIR/genome.fa
  READS_DIR="$GENOME_DIR/E${ERROR_RATE}_D${DEPTH}"
  READS=$READS_DIR/reads.fa
  ERROR_FLAGS="--p-mismatch $ERROR_RATE --p-gap-open $ERROR_RATE --p-gap-ext $ERROR_RATE --p-end 0.0"
  LOG="$READS_DIR/freqem_K${K}.tsv"

  $RELEASE -k $K benchmark $READS $GENOME -V 10 --init-state zero --parallel freq-em > $LOG
}

freqem     100 1 0 0.01 10  16
freqem     100 1 0 0.01 20  16
freqem     100 1 0 0.01 50  16
freqem     100 1 0 0.01 100 16

  # freqem     100 1 0 0.01 10  32
  # freqem     100 1 0 0.01 20  32
  # freqem     100 1 0 0.01 50  32
  # freqem     100 1 0 0.01 100 32
