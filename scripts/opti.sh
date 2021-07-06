export RUST_LOG=info
./target/debug/dbgphmm sample data/d1.fa --length 150 --n-reads 30 --start-from-head > data/r1.fa
./target/debug/dbgphmm read-stat data/d1.fa data/r1.fa
