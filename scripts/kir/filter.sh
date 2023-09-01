#!/bin/bash
set -Ceuo pipefail

REF='KIR_Primary.fa'
CUTOFF=2000
# TR region [690k:860k]
L=690000
R=860000

function filter () {
  READS=$1
  BASENAME=$(basename $READS)
  KEY=${BASENAME%.fastq.gz}
  PAF="$KEY.paf"

  echo $KEY $PAF

  # 1
  [ ! -e $PAF ] && minimap2 -t 16 -x map-pb -o $PAF $REF $READS
  cat $PAF | gawk "\$6 == \"LRC_KIR_Primary\" && \$8 <= $R && \$9 >= $L && \$10>=$CUTOFF" | cut -f1 | sort | uniq > $KEY.C$CUTOFF.ids

  # 2
  # seqkit fq2fa $READS | seqkit seq -m 1000 > ${ROOT}.all.fa
  seqkit faidx /work/ryought/mhc/data/KIR/filter/$KEY.all.fa --infile-list $KEY.C$CUTOFF.ids > $KEY.KIR.C$CUTOFF.fa
}

find /data/HG002/PacBio_CCS_15kb_20kb_chemistry2/*.fastq.gz | while read READS
do
  echo $READS
  filter $READS
done

# merge
cat *.KIR.C$CUTOFF.fa > reads.C$CUTOFF.fa
