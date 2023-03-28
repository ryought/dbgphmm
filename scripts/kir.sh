#!/bin/bash
#
#
cargo run --release --bin draft -- -k 20 -C 20 -p 0.01 -L 20000 -U 10000 -N 10 -E 300 -P 2 -H 0.01 --dataset-json-filename kir.json --gfa-filename kir.gfa --paths-filename kir.paths --dbg-filename kir.dbg

cargo run --release --bin draft -- -k 20 -C 20 -p 0.001 -L 20000 -U 10000 -N 10 -E 300 -P 2 -H 0.01 --dataset-json-filename kir2.json --gfa-filename kir2.gfa --paths-filename kir2.paths --dbg-filename kir2.dbg

cargo run --release --bin draft -- -k 16 -C 20 -p 0.01 -L 20000 -U 10000 -N 10 -E 300 -P 2 -H 0.01 --dataset-json-filename kir3.json --gfa-filename kir3.gfa --paths-filename kir3.paths --dbg-filename kir3.dbg
