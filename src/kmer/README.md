# kmer

## kmer definitions

- `VecKmer`
    normal kmer struct, it stores all bases in `Vec<u8>`
- `TinyKmer<K>`
    lightweight kmer struct by encoding bases in `u64`
- `Segment`
    Variable length k-mer struct, with `Vec<u8>`

## traits

- `KmerLike`
    A trait for fixed-length k-mer
- `KmerBase`
    A general trait for both variable-length and fixed-length k-mer


## Todos

- create variable k-mer struct
- add conversion tests
