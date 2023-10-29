#!/bin/bash
#
# Simulation functions
#

# to disable openblas threading
export OMP_NUM_THREADS=1

GEPARD_JAR=/home/ryought/data06/tools/gepard/gepard_fixed-20190402.jar
GEPARD_MAT=/home/ryought/data06/tools/gepard/edna.mat
# GEPARD_JAR=/Users/ryought/gepard/gepard_fixed-20190402.jar
# GEPARD_MAT=/Users/ryought/gepard/gepard_src/src/matrices/edna.mat

function gfa2fa () {
  # remove n gaps
  GFA=$1
  FA=${GFA/.gfa/.fa}
  awk '/^S/{print ">"$2;print $3}' $GFA | seqkit seq -g -G n | seqkit seq -m 1 > $FA
}

function map_to_genome() {
  GENOME=$1
  ASM=$2
  # minimap2 --secondary=no -c --cs -t4 -x asm20 $GENOME $ASM
  minimap2 -c --cs -t4 -x asm20 $GENOME $ASM
}

function genome_self_vs_self () {
  GENOME=$1
  # self-vs-self alignment
  map_to_genome $GENOME $GENOME > $GENOME.paf
  # awk '$1 != $6'
}

function summary_paf () {
  PAF=$1
  awk 'BEGIN{OFS="\t"} { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$13,$24 }' $1
}

function generate_svg () {
  GFA=$1
  PAF=$2
  python scripts/kir/graph_compare.py --min_identity 1.0 --spring_layout --hide_text $GFA $PAF
}

function generate_svg_for_dbgphmm () {
  GFA=$1
  PAF=$2
  EULER=$3
  python scripts/kir/graph_compare.py --min_identity 1.0 --hide_text --euler_fa $EULER $GFA $PAF
}

function gepard () {
  java -cp $GEPARD_JAR org.gepard.client.cmdline.CommandLine \
    -seq1 $1 -seq2 $2 -matrix $GEPARD_MAT -outfile $3 -maxwidth 1000 -maxheight 1000
}

function asm_eval () {
  GENOME=$1
  ASM=$2
  GFA=$3

  map_to_genome $GENOME $ASM > $ASM.paf

  generate_svg $GFA $ASM.paf > $ASM.svg

  gepard $ASM $ASM $ASM.png
}

function run_hifiasm () {
  KEY=$1
  READ="$KEY/data.reads.fa"
  GENOME="$KEY/data.genome.fa"
  DIR="$KEY/hifiasm"
  mkdir -p $DIR

  hifiasm -o $DIR/out -t4 -f0 -i $READ 2> $DIR/log

  gfa2fa $DIR/out.bp.p_ctg.gfa
  gfa2fa $DIR/out.bp.p_utg.gfa
  gfa2fa $DIR/out.bp.hap1.p_ctg.gfa
  gfa2fa $DIR/out.bp.hap2.p_ctg.gfa
  cat $DIR/out.bp.hap1.p_ctg.fa $DIR/out.bp.hap2.p_ctg.fa > $DIR/out.fa

  map_to_genome $GENOME $DIR/out.fa > $DIR/out.paf
  map_to_genome $GENOME $DIR/out.bp.p_utg.fa > $DIR/out.p_utg.paf

  gepard $DIR/out.fa $DIR/out.fa $DIR/out.fa.png
  gepard $GENOME $DIR/out.fa $DIR/out.fa.genome.png

  generate_svg $DIR/out.bp.p_utg.gfa $DIR/out.p_utg.paf > $DIR/out.svg
}

function run_verkko () {
  KEY=$1
  READ="$KEY/data.reads.fa"
  GENOME="$KEY/data.genome.fa"
  DIR="$KEY/verkko"
  mkdir -p $DIR

  verkko -d $DIR --hifi $READ

  if [ -e "$DIR/assembly.fasta" ]
  then
    asm_eval $GENOME $DIR/assembly.fasta $DIR/assembly.homopolymer-compressed.gfa
  fi
}

function run_lja () {
  KEY=$1
  READ="$KEY/data.reads.fa"
  GENOME="$KEY/data.genome.fa"
  DIR="$KEY/lja"
  mkdir -p $DIR

  lja -o $DIR --diploid --reads $READ

  map_to_genome $GENOME $DIR/assembly.fasta > $DIR/out.paf
  gepard $DIR/assembly.fasta $DIR/assembly.fasta $DIR/assembly.fasta.png
  gepard $GENOME $DIR/assembly.fasta $DIR/assembly.fasta.genome.png

  generate_svg $DIR/mdbg.gfa $DIR/out.paf > $DIR/out.svg
}

function evaluate_dbgphmm () {
  KEY=$1
  INFER_KEY=$2
  GENOME="$KEY/data.genome.fa"
  DIR="$KEY/dbgphmm"
  gfa2fa $DIR/$INFER_KEY.final.gfa
  map_to_genome $GENOME $DIR/$INFER_KEY.final.fa > $DIR/$INFER_KEY.final.paf

  generate_svg $DIR/$INFER_KEY.final.gfa $DIR/$INFER_KEY.final.paf > $DIR/$INFER_KEY.final.svg
  generate_svg_for_dbgphmm $DIR/$INFER_KEY.final.gfa $DIR/$INFER_KEY.final.paf $DIR/$INFER_KEY.final.euler.fa > $DIR/$INFER_KEY.final.euler.svg

  gepard $GENOME $DIR/$INFER_KEY.final.fa $DIR/$INFER_KEY.final.fa.genome.png
  gepard $GENOME $DIR/$INFER_KEY.final.euler.fa $DIR/$INFER_KEY.final.euler.fa.genome.png
}

function run_dbgphmm () {
  KEY=$1
  p=$2
  pi=$2
  DIR="$KEY/dbgphmm"
  pz=0.99
  INFER_KEY=pz${pz}_pi${pi}
  mkdir -p $DIR
  ./target/release/infer -t 32 -M 4 -k 40 -K 10000 -p $pi -e $p -s 5000 -I 50 --p0 $pz --dataset-json $KEY/data.json --output-prefix $DIR/$INFER_KEY

  evaluate_dbgphmm $KEY $INFER_KEY
}

function qsub_run_dbgphmm () {
  KEY=$1
  p=$2

  # execute run_dbgphmm $KEY $p in qsub
  # Key: replace / to -
  # -l hostname=z02 \
  qsub \
    -S /bin/bash -cwd -q all.q -j y -pe smp 32 \
    -N ${KEY//\//-} \
    << EOS
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/ryought/.miniconda3/lib:$LD_LIBRARY_PATH"
. scripts/sim.sh
echo "running $KEY with $p"
run_dbgphmm $KEY $p
EOS
}

function run_n4 () {
  p=0.0003

  for H in 0.01 0.001 0.0001
  do
    for H0 in 0.0002 0.0001
    do
      KEY="sim/n4_p${p}/H${H}_H0${H0}"
      mkdir -p $KEY
      echo $KEY

      # create dataset
      ./target/release/draft -k 40 -C 10 -L 10000 -p $p -M 4 -U 10000 -N 4 -E 2000 -H $H --H0 $H0 -P 2 --output-prefix $KEY/data --dataset-only

      genome_self_vs_self $KEY/data.genome.fa
      run_hifiasm $KEY
      run_lja $KEY

      # TODO activate miniconda to use verkko
      run_verkko $KEY

      # run_dbgphmm $KEY $p   # run locally
      qsub_run_dbgphmm $KEY $p   # run on cluster
      # evaluate_dbgphmm $KEY "pz0.99_pi0.0003"
    done
  done
}

function run_n10 () {
  p=0.0003

  for H in 0.01 0.001 0.0001
  do
    for H0 in 0.0002 0.0001
    do
      KEY="sim/n10_p${p}/H${H}_H0${H0}"
      mkdir -p $KEY
      echo $KEY

      # create dataset
      ./target/release/draft -k 40 -C 10 -L 10000 -p $p -M 4 -U 2000 -N 10 -E 2000 -H $H --H0 $H0 -P 2 --output-prefix $KEY/data --dataset-only

      # genome_self_vs_self $KEY/data.genome.fa
      run_hifiasm $KEY
      run_lja $KEY

      # TODO activate miniconda to use verkko
      # run_verkko $KEY

      qsub_run_dbgphmm $KEY $p   # run on cluster
      # run_dbgphmm $KEY $p   # run locally
      # evaluate_dbgphmm $KEY "pz0.99_pi0.0003"  # evaluate only
    done
  done
}

function manual_svg () {
  N=$1
  H=$2
  H0=$3
  HAPS0=$4
  HAPS1=$5
  KEY="sim/n${N}_p0.0003/H${H}_H0${H0}"
  PREFIX="$KEY/dbgphmm/pz0.99_pi0.0003"
  if [[ -n $HAPS0 && -n $HAPS1 ]]
  then
    python scripts/kir/graph_compare.py $PREFIX.final.gfa $PREFIX.final.paf --min_identity 1.0 --hide_text --haps $HAPS0 $HAPS1 > $PREFIX.final.svg
  else
    python scripts/kir/graph_compare.py $PREFIX.final.gfa $PREFIX.final.paf --min_identity 1.0 --hide_text --spring_layout > $PREFIX.final.svg
  fi
}

function svg_dbgphmm () {
  manual_svg 4 0.01   0.0001 4,1,2,5,2,3 4,6
  manual_svg 4 0.001  0.0001 1,8,4,0,6,3,4,5 1,9
  manual_svg 4 0.0001 0.0001
  manual_svg 4 0.01   0.0002 11,10,1,3,4 11,12
  manual_svg 4 0.001  0.0002 1,2,5,3 1,6
  manual_svg 4 0.0001 0.0002 0,8,5,6,10,7,9,2,14,13 0,12,5,15,10,7,9,2,1
}

function svg_hifiasm () {
  p=0.0003

  for n in 4 10
  do
    for H in 0.01 0.001 0.0001
    do
      for H0 in 0.0002 0.0001
      do
        KEY="sim/n${n}_p${p}/H${H}_H0${H0}"
        echo $KEY
        DIR="$KEY/hifiasm"
        python scripts/kir/graph_compare.py --min_identity 1.0 --spring_layout $DIR/out.bp.p_utg.gfa $DIR/out.p_utg.paf --hide_text > $DIR/out.svg
      done
    done
  done
}
