#!/bin/bash
export OMP_NUM_THREADS=1

G=50000

function gfa2fa () {
  GFA=$1
  FA=${GFA/.gfa/.fa}
  awk '/^S/{print ">"$2;print $3}' $GFA > $FA
}

function summary_paf () {
  PAF=$1
  awk 'BEGIN{OFS="\t"} { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$13,$24 }' $1
}

function run_hifiasm () {
  KEY=$1
  READ="t/$KEY/data.reads.fa"
  GENOME="t/$KEY/data.genome.fa"
  DIR="t/$KEY/hifiasm_opt"
  mkdir -p $DIR
  # hifiasm -o $DIR/out -t4 -f0 $READ 2> $DIR/log
  hifiasm -o $DIR/out -t4 -f0 --hg-size 70k -D 50 -i $READ 2> $DIR/log

  gfa2fa $DIR/out.bp.p_ctg.gfa
  gfa2fa $DIR/out.bp.p_utg.gfa
  gfa2fa $DIR/out.bp.hap1.p_ctg.gfa
  gfa2fa $DIR/out.bp.hap2.p_ctg.gfa
  cat $DIR/out.bp.hap1.p_ctg.fa $DIR/out.bp.hap2.p_ctg.fa > $DIR/out.fa

  minimap2 --secondary=no -c --cs -t4 -x asm20 $GENOME $DIR/out.fa > $DIR/out.paf
  minimap2 --secondary=no -c --cs -t4 -x asm20 $GENOME $DIR/out.bp.p_utg.fa > $DIR/out.p_utg.paf
}



function run_v0 () {
  for C in 3 5 10
  do
    for U in 10000 5000 2500 1000 500
    do
      N=$(( $G / $U ))
      for H in 0.02 0.01 0.001
      do
        KEY="U${U}_N${N}_H${H}_C${C}"
        echo "$C $U $N $H $KEY"
        # mkdir -p t/$KEY

        # create dataset
        ./target/release/draft -k 40 -C $C -L 10000 -p 0.001 -m 1 -M 2 -U $U -N $N -E 10000 -H $H --H0 $H -P 2 --output-prefix t/$KEY/data --dataset-only

        # run hifiasm
        # run_hifiasm $KEY

        summary_paf "t/$KEY/hifiasm/out.p_utg.paf"
        # summary_paf "t/$KEY/hifiasm_opt/out.paf"
      done
    done
  done
  # ./target/release/infer -k 40 -p 0.00001 -K 10000 -e 0.001 -I 50 -s 10000 --dbg t/t8.dbg --dataset-json t/t8.json --output-prefix t/t8
}


function run_dbgphmm () {
  C=10
  U=5000
  H=0.001
  N=$(( $G / $U ))
  KEY="U${U}_N${N}_H${H}_C${C}"
  ./target/release/draft -k 40 -C $C -L 10000 -p 0.001 -M 4 -U $U -N $N -E 10000 -H $H --H0 $H -P 2 --output-prefix t/$KEY/data
  ./target/release/infer --dbg t/$KEY/data.dbg -K 10000 -p 0.00001 -e 0.001 -s 10000 -I 50 --dataset-json t/$KEY/data.json --output-prefix t/$KEY/v0
}
