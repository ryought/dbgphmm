export OMP_NUM_THREADS=1
./target/release/draft -k 40 --genome-fasta kir/haps01.fa -C 20 -L 20000 -p 0.001 -U 10 -N 10 -E 10 -H 0.0 --H0 0.0 -P 2 --output-prefix kir/haps/v0
./target/release/infer -k 40 -K 20000 -p 0.00001 -e 0.001 -I 50 -s 10000 --dataset-json kir/haps/v0.json --output-prefix kir/haps/v0
