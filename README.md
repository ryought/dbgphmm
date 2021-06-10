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
