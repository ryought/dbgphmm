#!/bin/bash
# fetch T2T HG002 assembly from AWS and determine corresponding region of LRC_KIR_Primary:680000-860000

# Input:
# Fasta file including LRC_KIR_Primary:680000-860000 (G0 haplotype of KIR)
HAPS='haps01.fa'

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/assemblies/drafts/assembly.v0.7.fasta
seqkit faidx assembly.v0.7.fasta chr19_MATERNAL chr19_PATERNAL > T2T_chr19.fa

# 1000bp prefix/suffix of hg38 KIR region
seqkit faidx $HAPS 'LRC_KIR_Primary:680000-860000':1-1000 > prefix.fa
seqkit faidx $HAPS 'LRC_KIR_Primary:680000-860000':-1000--1 > suffix.fa

minimap2 -c --cs -t 16 -x asm20 T2T_chr19.fa prefix.fa > prefix.paf
minimap2 -c --cs -t 16 -x asm20 T2T_chr19.fa suffix.fa > suffix.paf

cut -f6-9 prefix.paf
# chr19_MATERNAL  61325064        57332152        57333152
# chr19_PATERNAL  61370595        57430156        57431156

cut -f6-9 suffix.paf
# chr19_MATERNAL  61325064        57586627        57587628
# chr19_PATERNAL  61370595        57636207        57637208

MAT_L=57332152
MAT_R=57587628
PAT_L=57430156
PAT_R=57637208

seqkit faidx assembly.v0.7.fasta \
  chr19_MATERNAL:$(($MAT_L + 1))-$MAT_R \
  chr19_PATERNAL:$(($PAT_L + 1))-$PAT_R \
  > KIR.fa
