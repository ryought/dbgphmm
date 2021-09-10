#!/bin/bash
# count-base
# good?

# dbgphmm -k 8 generate -l 300 > g1.fa
dbgphmm -k 8 stat g1.fa > g1.json
dbgphmm -k 400 sample g1.fa --length 1000 --n-reads 10 --start-from-head > g1_r1.fa
# dbgphmm -k 8 --p-mismatch 0 --p-gap-open 0 --p-gap-ext 0 --p-end 0.0 sample g1.fa --length 10000 --n-reads 10 --start-from-head > g1_r2.fa

# dbgphmm -k 8 benchmark g1_r1.fa g1.fa -V 10 --init-state zero freq-em
