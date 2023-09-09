#!/bin/bash
set -Ceuo pipefail

# Test 1: use
# REF='KIR_Primary.fa'
# at least 2000 bp
# CUTOFF=2000
# TR region
# L=680000
# R=860000


# Test 2
REF='chr19.fa'
CUTOFF=2000
# LRC_KIR region
START=54025634
END=55084318
# tandem repeat region in chr19
# 680 kb : 860 kb
L=$(( $START + 680000 ))
R=$(( $START + 860000 ))
# 54,705 kb : 54,905 kb
# 54,896 : 54,900 kb
# --start 54_705_634 --end 54_885_634
# 868 874


cd ~/work/dbgphmm/kir

function map () {
  READS=$1
  BASENAME=$(basename $READS)
  KEY=${BASENAME%.fastq.gz}
  PAF="$KEY.paf"

  echo $KEY $PAF
  # [ ! -e $PAF ] && minimap2 -t 16 -x map-pb -o $PAF $REF $READS
  minimap2 -t 16 -x asm20 -o $PAF $REF $READS
}

function filter () {
  READS=$1
  BASENAME=$(basename $READS)
  KEY=${BASENAME%.fastq.gz}
  PAF="$KEY.paf"

  # forward
  echo forward
  gawk "\$8 <= $R && \$9 >= $L && \$10>=$CUTOFF && \$5 == \"+\"" $PAF > $KEY.f.paf
  cut -f1 $KEY.f.paf | sort | uniq > $KEY.f.ids
  # backward
  echo backward
  gawk "\$8 <= $R && \$9 >= $L && \$10>=$CUTOFF && \$5 == \"-\"" $PAF > $KEY.b.paf
  cut -f1 $KEY.b.paf | sort | uniq > $KEY.b.ids


  # 2
  # seqkit fq2fa $READS | seqkit seq -m 1000 > ${ROOT}.all.fa
  # forward
  seqkit faidx /work/ryought/mhc/data/KIR/filter/$KEY.all.fa --infile-list $KEY.f.ids > $KEY.f.fa
  seqkit faidx /work/ryought/mhc/data/KIR/filter/$KEY.all.fa --infile-list $KEY.b.ids | seqkit seq -rpv -t dna > $KEY.b.fa
}

function mapping_with_igv () {
  READS=$1
  echo 'mapping'
  minimap2 -c -x asm20 --cs -a $REF $READS > $READS.sam
  echo 'converting bam'
  samtools view -bS $READS.sam | samtools sort > $READS.bam
  samtools index $READS.bam
}

function merge () {
  cat *.f.fa *.b.fa > merged.fa
}

find /data/HG002/PacBio_CCS_15kb_20kb_chemistry2/*.fastq.gz | while read READS
do
  echo $READS
  # map $READS
  # filter $READS
done
# merge
# mapping_with_igv merged.fa

python scripts/kir/sam_to_fa.py merged.fa.sam merged.fa --start 54705634 --end 54885634 > trimed.fa
mapping_with_igv trimed.fa

# subsample
seqkit sample -p 0.5 -s 0 trimed.fa > 20x.fa
seqkit sample -p 0.25 -s 1 trimed.fa > 10x.fa
