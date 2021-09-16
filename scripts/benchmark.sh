#$ -S /bin/bash
#$ -cwd
#$ -V
#$ -o ./out
#$ -j y
#$ -pe smp 6

# to submit
# $ qsub -v L=50 -v S=0 -N "dbgphmm_L50_S0" scripts/benchmark.sh
set -u

DEBUG='cargo run --'
RELEASE='./target/release/dbgphmm'

function run () {
  LENGTH=$1
  N_HAP=$2
  SEED=$3
  ERROR_RATE=$4
  DEPTH=$5
  K=$6
  FLAGS=$7
  LOG_KEY=$8

  GENOME_DIR="data/L${LENGTH}_N${N_HAP}_S${SEED}"
  GENOME=$GENOME_DIR/genome.fa
  READS_DIR="$GENOME_DIR/E${ERROR_RATE}_D${DEPTH}"
  READS=$READS_DIR/reads.fa
  ERROR_FLAGS="--p-mismatch $ERROR_RATE --p-gap-open $ERROR_RATE --p-gap-ext $ERROR_RATE --p-end 0.0"
  LOG="$READS_DIR/${LOG_KEY}_K${K}.tsv"
  ERR="$READS_DIR/${LOG_KEY}_K${K}.err"

  $RELEASE -k $K benchmark $READS $GENOME -V 10 --parallel $FLAGS 1> $LOG 2> $ERR
}

function main () {
  L=$1
  S=$2
  for E in 0.01 0.02 0.05
  do
    for D in 10 50 100
    do
      for K in 8 16 24 32
      do
        run $L 1 $S $E $D $K '--init-state zero freq-em' 'freqem2'
        run $L 1 $S $E $D $K '--init-state read-count full-em --depth-scheduler linear-gradient -I 20' 'fullemlingrad2'
        run $L 1 $S $E $D $K '--init-state read-count full-em --depth-scheduler constant -I 20' 'fullemconst2'
        # run $L 1 $S $E $D $K '--init-state zero grad -I 10' 'grad_f0_I10'
        # run $L 1 $S $E $D $K '--init-state true grad -I 10' 'grad_ft_I10'
      done
    done
  done
}

main $L $S
