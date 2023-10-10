#!/bin/bash
#$ -S /bin/bash
#$ -N dbgphmm
#$ -cwd
#$ -q all.q
#$ -j y
#$ -l hostname=z02
##$ -pe smp 128

# cd /work/00/gg57/j29006/dbgphmm
# module load python/3.7.3
export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH="/home/ryought/.miniconda3/lib:$LD_LIBRARY_PATH"

. scripts/compare.sh
run_dbgphmm_u2k
