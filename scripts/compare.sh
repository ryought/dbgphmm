#!/bin/bash
export OMP_NUM_THREADS=1

G=50000

function gfa2fa () {
  GFA=$1
  FA=${GFA/.gfa/.fa}
  awk '/^S/{print ">"$2;print $3}' $GFA > $FA
}

function map_to_genome() {
  GENOME=$1
  ASM=$2
  # minimap2 --secondary=no -c --cs -t4 -x asm20 $GENOME $ASM
  minimap2 -c --cs -t4 -x asm20 $GENOME $ASM
}

function summary_paf () {
  PAF=$1
  awk 'BEGIN{OFS="\t"} { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$13,$24 }' $1
}

function generate_svg () {
  GFA=$1
  PAF=$2
  # --draw_mismatch_primary_only
  python scripts/kir/graph_compare.py --draw_mismatch_threshold 100 $GFA $PAF
}

function generate_svg_primary () {
  GFA=$1
  PAF=$2
  python scripts/kir/graph_compare.py --draw_mismatch_primary_only $GFA $PAF
}

function gepard () {
  java -cp /home/ryought/data06/tools/gepard/gepard_fixed-20190402.jar org.gepard.client.cmdline.CommandLine \
    -seq1 $1 -seq2 $2 -matrix /home/ryought/data06/tools/gepard/edna.mat -outfile $3 -maxwidth 1000 -maxheight 1000
}

function run_hifiasm () {
  KEY=$1
  READ="t/$KEY/data.reads.fa"
  GENOME="t/$KEY/data.genome.fa"
  DIR="t/$KEY/hifiasm"
  # DIR="t/$KEY/hifiasm_opt"
  mkdir -p $DIR
  hifiasm -o $DIR/out -t4 -f0 -i $READ 2> $DIR/log
  # hifiasm -o $DIR/out -t4 -f0 --hg-size 70k -D 50 -i $READ 2> $DIR/log

  gfa2fa $DIR/out.bp.p_ctg.gfa
  gfa2fa $DIR/out.bp.p_utg.gfa
  gfa2fa $DIR/out.bp.hap1.p_ctg.gfa
  gfa2fa $DIR/out.bp.hap2.p_ctg.gfa
  cat $DIR/out.bp.hap1.p_ctg.fa $DIR/out.bp.hap2.p_ctg.fa > $DIR/out.fa

  map_to_genome $GENOME $DIR/out.fa > $DIR/out.paf
  map_to_genome $GENOME $DIR/out.bp.p_utg.fa > $DIR/out.p_utg.paf

  generate_svg $DIR/out.bp.p_utg.gfa $DIR/out.p_utg.paf > $DIR/out.svg
  generate_svg_primary $DIR/out.bp.p_utg.gfa $DIR/out.p_utg.paf > $DIR/out.primary.svg

  gepard $DIR/out.fa $DIR/out.fa $DIR/out.fa.png
}

function run_verkko () {
  KEY=$1
  READ="t/$KEY/data.reads.fa"
  GENOME="t/$KEY/data.genome.fa"
  DIR="t/$KEY/verkko"
  mkdir -p $DIR

  verkko -d $DIR --hifi $READ

  map_to_genome $GENOME $DIR/assembly.fasta > $DIR/out.paf
  # FIXME
  generate_svg $DIR/assembly.homopolymer-compressed.gfa $DIR/out.paf > $DIR/out.svg
}

function run_lja () {
  KEY=$1
  READ="t/$KEY/data.reads.fa"
  GENOME="t/$KEY/data.genome.fa"
  DIR="t/$KEY/lja"
  mkdir -p $DIR

  lja -o $DIR --diploid --reads $READ

  map_to_genome $GENOME $DIR/assembly.fasta > $DIR/out.paf
  generate_svg $DIR/mdbg.gfa $DIR/out.paf > $DIR/out.svg
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

function run_v1 () {
  C=10
  for N in 2 5 10
  do
    U=$(( $G / $N ))
    for H in 0.002 0.001
    do
      for H0 in 0.001 0
      do
        KEY="U${U}_N${N}_H${H}_H0${H0}_C${C}"
        echo "$C $U $N $H $KEY"
        mkdir -p t/$KEY

        # create dataset
        ./target/release/draft -k 40 -C $C -L 10000 -p 0.001 -m 1 -M 2 -U $U -N $N -E 10000 -H $H --H0 $H0 -P 2 --output-prefix t/$KEY/data --dataset-only

        # run hifiasm
        run_hifiasm $KEY

        # summary_paf "t/$KEY/hifiasm/out.p_utg.paf"
        # summary_paf "t/$KEY/hifiasm_opt/out.paf"
      done
    done
  done
}

function run_v2 () {
  C=10
  p=0.0003
  # p=0.001
  U=2000
  L=10000
  # for N in 2 3 4 5
  # for N in 2 3 4
  for N in 10
  do
    for H in 0.01 0.001 0.0001
    do
      for H0 in 0.001 0.0001
      do
        KEY="U${U}_N${N}_H${H}_H0${H0}_C${C}_p${p}"
        mkdir -p t/$KEY
        echo $KEY

        # create dataset
        ./target/release/draft -k 40 -C $C -L $L -p $p -m 1 -M 2 -U $U -N $N -E 2000 -H $H --H0 $H0 -P 2 --output-prefix t/$KEY/data --dataset-only
        # ./target/release/draft -k 40 -C $C -L 10000 -p $p -M 4 -U $U -N $N -E 10000 -H $H --H0 $H0 -P 2 --output-prefix t/$KEY/data

        # self-vs-self alignment
        # map_to_genome t/$KEY/data.genome.fa t/$KEY/data.genome.fa > t/$KEY/data.genome.paf
        # awk '$1 != $6' t/$KEY/data.genome.paf

        # run hifiasm
        run_hifiasm $KEY
      done
    done
  done
}


function run_dbgphmm () {
  C=10
  U=5000
  H=0.001
  N=$(( $G / $U ))
  KEY="U${U}_N${N}_H${H}_C${C}"
  ./target/release/draft -k 40 -C $C -L 10000 -p 0.001 -M 4 -U $U -N $N -E 10000 -H $H --H0 $H -P 2 --output-prefix t/$KEY/data
  ./target/release/infer --dbg t/$KEY/data.dbg -K 10000 -p 0.00001 -e 0.001 -s 10000 -I 50 --dataset-json t/$KEY/data.json --output-prefix t/$KEY/dbgphmm

  gfa2fa t/$KEY/dbgphmm/*.final.gfa
  map_to_genome $GENOME $DIR/out.fa > $DIR/out.paf
}

function benchmark_dbgphmm () {
  KEY=$1
  INFER_KEY=$2
  GENOME="t/$KEY/data.genome.fa"
  DIR="t/$KEY/dbgphmm"
  gfa2fa $DIR/$INFER_KEY.final.gfa
  map_to_genome $GENOME $DIR/$INFER_KEY.final.fa > $DIR/$INFER_KEY.final.paf
  generate_svg_primary $DIR/$INFER_KEY.final.gfa $DIR/$INFER_KEY.final.paf > $DIR/$INFER_KEY.final.svg
}

# run_v2
# run_v1
# run_hifiasm "U10000_N5_H0.001_H00.0_C10_p0.0003"
# run_lja "U10000_N5_H0.001_H00.0_C10_LJA"
# run_hifiasm "U10000_N3_H0.001_H00.0_C10_p0.001"
# run_hifiasm "U25000_N2_H0.001_H00.0_C10_p0.001"
# run_hifiasm "U25000_N2_H0.001_H00.001_C10"
