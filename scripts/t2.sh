#!/bin/bash
#$ -S /bin/bash
#$ -N dbgphmm
#$ -cwd
#$ -q all.q
#$ -j y
#$ -l hostname=g01
##$ -pe smp 128

# cd /work/00/gg57/j29006/dbgphmm
# module load python/3.7.3
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/ryought/.miniconda3/lib:$LD_LIBRARY_PATH"

C=10
# 0.01%
H=0.0001
# H0=0.0
H0=0.0001
# G=50000
N=3
# U=$(( $G / $N ))
U=10000
# 0.03% for LJA
# p=0.0003
# 0.1% (~hifi)
p=0.001
KEY="U${U}_N${N}_H${H}_H0${H0}_C${C}_p${p}"

# U10000_N3_H0.0001_H00.0002_C10_p0.001

# 0.1%
# pi=0.0001
# pz=0.9999
pi=0.001
pz=0.99

echo running $KEY

mkdir -p t/$KEY/dbgphmm
./target/release/draft -k 40 -C $C -L 10000 -p $p -M 4 -U $U -N $N \
  -E 10000 -H $H --H0 $H0 -P 2 --output-prefix t/$KEY/data
  # --read-seed 5 \
./target/release/infer --dbg t/$KEY/data.dbg -K 10000 -p $pi -e $p --p0 $pz -s 10000 -I 50 --dataset-json t/$KEY/data.json --output-prefix t/$KEY/dbgphmm/pz${pz}_pi${pi}
