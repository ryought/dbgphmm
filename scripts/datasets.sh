#!/bin/bash
set -u

DEBUG='cargo run --'
RELEASE='./target/release/dbgphmm'

function create_genome () {
  LENGTH=$1
  N_HAP=$2
  SEED=$3

  # genome
  GENOME_DIR="data/L${LENGTH}_N${N_HAP}_S${SEED}"
  mkdir -p $GENOME_DIR
  GENOME=$GENOME_DIR/genome.fa
  $RELEASE generate --length $LENGTH --seed $SEED > $GENOME
  # TODO generate haplotype by inserting random mutation
  # TODO more complex genomes, such as duplication etc
}

function create_reads () {
  LENGTH=$1
  N_HAP=$2
  SEED=$3
  ERROR_RATE=$4
  DEPTH=$5

  GENOME_DIR="data/L${LENGTH}_N${N_HAP}_S${SEED}"
  GENOME=$GENOME_DIR/genome.fa

  # reads
  READS_DIR="$GENOME_DIR/E${ERROR_RATE}_D${DEPTH}"
  mkdir -p $READS_DIR
  READS=$READS_DIR/reads.fa
  ERROR_FLAGS="--p-mismatch $ERROR_RATE --p-gap-open $ERROR_RATE --p-gap-ext $ERROR_RATE --p-end 0.0"
  # XXX length is set to large enough to span whole genome
  $RELEASE -k 32 $ERROR_FLAGS sample $GENOME --length 10000 --n-reads $DEPTH --start-from-head > $READS
}

function main () {
  for L in 50 100 500
  do
    for S in 0 1 2
    do
      create_genome $L 1 $S

      for E in 0.01 0.02 0.05
      do
        for D in 10 50 100
        do
          echo $L $S $E $D
          create_reads  $L 1 $S $E $D
        done
      done
    done
  done
}

main
