#!/bin/sh
#PJM -g gg57
#PJM -L rscgrp=regular
#PJM -L node=1
#PJM -L elapse=48:00:00
#PJM -j

# pjsub -N lja -j scripts/kir/lja.sh

cd /work/00/gg57/j29006/dbgphmm

module unload gcc/4.8.5 intel
module load cmake/3.14.5 gcc/7.5.0

# /work/gg57/share/ryought-tools/LJA/bin/lja -o kir/haps/LJA/10x --reads kir/haps/10x.reads.fa
# /work/gg57/share/ryought-tools/LJA/bin/lja -o kir/haps/LJA/20x --reads kir/haps/20x.reads.fa --diploid

# /work/gg57/share/ryought-tools/LJA/bin/lja -o kir/hifi/LJA/20x --reads kir/hifi/20x.fa --diploid
/work/gg57/share/ryought-tools/LJA/bin/lja -o kir/hifi/LJA/10x --reads kir/hifi/10x.fa --diploid
