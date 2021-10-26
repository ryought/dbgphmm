cd data/rep

dbgphmm -k 8 generate -l 50 -s 1 > prefix.fa
dbgphmm -k 8 generate -l 50 -s 2 > suffix.fa
dbgphmm -k 8 generate -l 50 -s 3 > unit01.fa
seqkit concat -w 10000 prefix.fa unit01.fa unit01.fa unit01.fa suffix.fa 2> /dev/null > rep01.fa
seqkit concat -w 10000 prefix.fa unit01.fa unit02.fa unit03.fa suffix.fa 2> /dev/null > rep02.fa
