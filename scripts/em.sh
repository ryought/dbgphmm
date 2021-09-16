#!/bin/bash
# full-em test

# genome
dbgphmm -k 8 generate -l 10 > g0.fa
dbgphmm -k 8 generate -l 11 > g4.fa
dbgphmm -k 8 generate -l 40 > g1.fa
dbgphmm -k 8 stat g1.fa > g1.json

seqkit concat -w 10000 g1.fa g1.fa 2> /dev/null > g2.fa
seqkit concat -w 10000 g0.fa g1.fa g1.fa g4.fa 2> /dev/null > g3.fa
# dbgphmm -k 60 stat g2.fa > g2.json
dbgphmm -k 8 stat g2.fa > g2.json
dbgphmm -k 8 stat g3.fa > g3.json

# read
# dbgphmm -k 30 sample g1.fa --length 1000 --n-reads 10 --start-from-head > g1_r1.fa
# dbgphmm -k 60 sample g2.fa --length 1000 --n-reads 10 --start-from-head > g2_r1.fa
dbgphmm -k 60 sample g3.fa --length 1000 --n-reads 10 --start-from-head > g3_r1.fa

# error-free read
# dbgphmm -k 8 --p-mismatch 0 --p-gap-open 0 --p-gap-ext 0 --p-end 0.0 sample g1.fa --length 10000 --n-reads 10 --start-from-head > g1_r2.fa

# dbgphmm -k 8 benchmark g1_r1.fa g1.fa -V 10 --init-state read-count full-em --depth-scheduler linear-gradient
# dbgphmm -k 8 benchmark g1_r1.fa g1.fa -V 10 --init-state zero freq-em
# dbgphmm -k 8 benchmark g2_r1.fa g2.fa -V 10 --init-state read-count full-em --depth-scheduler linear-gradient

dbgphmm -k 8 benchmark g3_r1.fa g3.fa --parallel -V 10 --init-state read-count full-em --depth-scheduler linear-gradient
# dbgphmm -k 8 benchmark g3_r1.fa g3.fa -V 10 --init-state zero freq-em
