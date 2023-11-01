# dbgphmm

Bayesian genome assembler using de Bruijn graph and profile HMM

## Build

Requires nightly rust

```
cargo build --release
./target/release/dbgphmm -h
```

## Simulation experiments

```
./target/release/draft -h
./target/release/infer -h
```

## Outputs

* `${prefix}.k${k}.gfa`
* `${prefix}.k${k}.inspect`
* `${prefix}.k${k}.dbg`
* `${prefix}.k${k}.post`

* `${prefix}.final.gfa`
* `${prefix}.final.euler.fa`


### INSPECT

### DBG



## Citation

Bayesian genome assembly of segmental duplications by inferring k-mer copy numbers in de Bruijn graphs


## WIP: python binding

```
source .env/bin/activate
maturin build
pip install --force-reinstall target/wheels/dbgphmm-0.1.0-cp310-cp310-macosx_11_0_arm64.whl
python -c 'import dbgphmm; print(repr(dbgphmm.sum_as_string(1, 2)));'
```
