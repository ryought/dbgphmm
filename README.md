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
