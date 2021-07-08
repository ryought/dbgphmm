
./target/release/dbgphmm -k 8 optimize data/r1.fa --true-dbg-fa data/d1.fa -M 10 -V 10 -I 100 -R 0.8 --parallel --start-from-true-copy-nums > data/r1.tsv
./target/release/dbgphmm -k 8 stat data/r1.fa > data/r1.json
