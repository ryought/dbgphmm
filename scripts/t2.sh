#!/bin/bash
#$ -S /bin/bash
#$ -N U2500
#$ -cwd
#$ -q centos7.q
#$ -j y

# cd /work/00/gg57/j29006/dbgphmm
# module load python/3.7.3
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/ryought/.miniconda3/lib:$LD_LIBRARY_PATH"

C=10
H=0.001
G=50000
# U=1000
# U=10000
U=2500
N=$(( $G / $U ))
KEY="U${U}_N${N}_H${H}_C${C}"

echo running $KEY

mkdir -p t/$KEY
./target/release/draft -k 40 -C $C -L 10000 -p 0.001 -M 4 -U $U -N $N -E 10000 -H $H --H0 $H -P 2 --output-prefix t/$KEY/data
./target/release/infer --dbg t/$KEY/data.dbg -K 10000 -p 0.00001 -e 0.001 -s 10000 -I 50 --dataset-json t/$KEY/data.json --output-prefix t/$KEY/v0
