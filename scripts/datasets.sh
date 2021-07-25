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

# create_genome 100 1 0
  create_reads  100 1 0 0.00 10

  create_reads  100 1 0 0.01 10
  create_reads  100 1 0 0.01 20
  create_reads  100 1 0 0.01 50
  create_reads  100 1 0 0.01 100

  create_reads  100 1 0 0.02 10
# create_genome 100 1 1
