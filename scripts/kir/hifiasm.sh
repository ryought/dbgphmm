#!/bin/sh
#PJM -g gg57
#PJM -L rscgrp=regular
#PJM -L node=1
#PJM -L elapse=48:00:00
#PJM -j

# pjsub -N hifiasm -j scripts/kir/hifiasm.sh

cd /work/00/gg57/j29006/dbgphmm
# haps
# /work/gg57/share/ryought-tools/hifiasm/hifiasm -o kir/haps/hifiasm/10x -t32 -f0 kir/haps/10x.reads.fa
# hifi
# /work/gg57/share/ryought-tools/hifiasm/hifiasm -o kir/hifi/hifiasm/10x -t32 -f0 kir/hifi/10x.fa
/work/gg57/share/ryought-tools/hifiasm/hifiasm -o kir/hifi/hifiasm/20x -t32 -f0 kir/hifi/20x.fa
