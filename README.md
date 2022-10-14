# dbgphmm - de bruijn graph + profile HMM
package

## crates
- src/main.rs (binary)
- src/lib.rs
- src/bin/hoge.rs (binary)

`mod my;`: search for `my.rs` or `my/mod.rs` from the file directory
    - main.rs -> src/my.rs or src/my/mod.rs
    - hoge.rs -> src/hoge/my.rs
`mod my {}`

## python binding

```
source .env/bin/activate
maturin build
pip install --force-reinstall target/wheels/dbgphmm-0.1.0-cp310-cp310-macosx_11_0_arm64.whl
python -c 'import dbgphmm; print(repr(dbgphmm.sum_as_string(1, 2)));'
```

## log level
- DEBUG only for small data
- INFO dev
- WARN production

## development
```
# For debug
cargo build
./target/debug/dbgphmm
./target/debug/hoge
or
cargo run

# For release
cargo build --release
./target/release/dbgphmm
```

### run test with stdout
```sh
cargo test -- --nocapture
```

### set log level
```sh
RUST_LOG=info cargo run -- sandbox
RUST_LOG=dbgphmm::cli::info cargo run -- sandbox
```

### automatic build update
```sh
cargo watch --ignore scripts/ -x 'build --release'
```

### benchmark


## Datasets
```sh
data/L${LEN}_N${N_HAP}_S${SEED}/genome.fa

data/L${LEN}_N${N_HAP}_S${SEED}/E${ERROR_RATE}_D${DEPTH}/reads.fa

data/L${LEN}_N${N_HAP}_S${SEED}/E${ERROR_RATE}_D${DEPTH}/stats.json
data/L${LEN}_N${N_HAP}_S${SEED}/E${ERROR_RATE}_D${DEPTH}/grad.tsv
```
