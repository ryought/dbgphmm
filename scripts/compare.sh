#!/bin/bash
export OMP_NUM_THREADS=1

G=50000

. scripts/sim.sh


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
  # p=0.0003
  p=0.001
  # U=2000
  U=10000
  L=10000
  E=2000
  # for N in 2 3 4 5
  # for N in 10
  for N in 3
  do
    for H in 0.01 0.001 0.0001
    do
      for H0 in 0.0 0.0001
      do
        KEY="U${U}_N${N}_H${H}_H0${H0}_C${C}_p${p}_E${E}"
        mkdir -p t/$KEY
        echo $KEY

        # create dataset
        ./target/release/draft -k 40 -C $C -L $L -p $p -m 1 -M 2 -U $U -N $N -E $E -H $H --H0 $H0 -P 2 --output-prefix t/$KEY/data --dataset-only
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

function run_v2b () {
  C=10
  p=0.0003
  U=10000
  L=10000
  E=2000
  # for N in 3 4
  for N in 3 4
  do
    # for H in 0.0001 0.001 0.01
    # for H in 0.0001
    for H in 0.01
    do
      for H0 in 0.0002 0.0001
      do
        KEY="U${U}_N${N}_H${H}_H0${H0}_C${C}_p${p}"
        mkdir -p t/$KEY
        echo $KEY

        # create dataset
        ./target/release/draft -k 40 -C $C -L $L -p $p -m 1 -M 2 -U $U -N $N -E $E -H $H --H0 $H0 -P 2 --output-prefix t/$KEY/data --dataset-only
        # ./target/release/draft -k 40 -C $C -L 10000 -p $p -M 4 -U $U -N $N -E 10000 -H $H --H0 $H0 -P 2 --output-prefix t/$KEY/data

        # self-vs-self alignment
        map_to_genome t/$KEY/data.genome.fa t/$KEY/data.genome.fa > t/$KEY/data.genome.paf
        awk '$1 != $6' t/$KEY/data.genome.paf

        # run hifiasm
        run_hifiasm $KEY
        # DIR="t/$KEY/hifiasm"
        # gepard $DIR/out.fa $DIR/out.fa $DIR/out.fa.png
      done
    done
  done
}

function run_dbgphmm_u2k () {
  H=0.01
  H0=0.0001
  # p=0.0003
  # or
  p=0.001
  KEY="U2kN10_H${H}_H0${H0}_p${p}"
  mkdir -p t/$KEY
  # pi=0.0001
  pi=0.001
  pz=0.99
  echo run $KEY
  ./target/release/draft -k 40 -C 10 -L 10000 -p $p -M 4 -U 2000 -N 10 -E 2000 -H $H --H0 $H0 -P 2 --output-prefix t/$KEY/data
  mkdir -p t/$KEY/dbgphmm
  ./target/release/infer --dbg t/$KEY/data.dbg -K 10000 -p $pi -e $p -s 4000 -I 50 --p0 $pz --dataset-json t/$KEY/data.json --output-prefix t/$KEY/dbgphmm/pz${pz}_pi${pi}

  gfa2fa t/$KEY/dbgphmm/*.final.gfa
  map_to_genome $GENOME $DIR/$INFER_KEY.final.fa > $DIR/$INFER_KEY.final.paf
  generate_svg_primary $DIR/$INFER_KEY.final.gfa $DIR/$INFER_KEY.final.paf > $DIR/$INFER_KEY.final.primary.svg
  generate_svg $DIR/$INFER_KEY.final.gfa $DIR/$INFER_KEY.final.paf > $DIR/$INFER_KEY.final.svg
}

function run_dbgphmm_u10k () {
  H=0.01    # 1%
  H0=0.0002 # 0.01% or 0.02%
  p=0.0003  # 0.03%
  # or
  # p=0.001
  KEY="U10kN4_H${H}_H0${H0}_p${p}"
  mkdir -p t/$KEY
  pi=0.0003
  pz=0.99
  INFER_KEY=pz${pz}_pi${pi}
  echo run $KEY $INFER_KEY
  ./target/release/draft -k 40 -C 10 -L 10000 -p $p -M 4 -U 10000 -N 4 -E 2000 -H $H --H0 $H0 -P 2 --output-prefix t/$KEY/data
  mkdir -p t/$KEY/dbgphmm
  ./target/release/infer --dbg t/$KEY/data.dbg -K 10000 -p $pi -e $p -s 10000 -I 50 --p0 $pz --dataset-json t/$KEY/data.json --output-prefix t/$KEY/dbgphmm/$INFER_KEY

  gfa2fa t/$KEY/dbgphmm/*.final.gfa
  DIR=t/$KEY/dbgphmm
  map_to_genome $GENOME $DIR/$INFER_KEY.final.fa > $DIR/$INFER_KEY.final.paf
  generate_svg_primary $DIR/$INFER_KEY.final.gfa $DIR/$INFER_KEY.final.paf > $DIR/$INFER_KEY.final.primary.svg
  generate_svg $DIR/$INFER_KEY.final.gfa $DIR/$INFER_KEY.final.paf > $DIR/$INFER_KEY.final.svg
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
