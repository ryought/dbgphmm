#![feature(test)]
#![feature(min_const_generics)]
pub mod cli;
pub mod compressed_dbg;
pub mod copy_nums;
pub mod cycles;
pub mod dbg;
pub mod distribution;
pub mod graph;
pub mod hmm;
pub mod io;
pub mod kmer;
pub mod mocks;
pub mod optimizer;
pub mod prob;
pub mod random_seq;
pub mod sparse;
pub mod stats;
pub mod veclike;

extern crate jemallocator;
#[global_allocator]
static ALLOC: jemallocator::Jemalloc = jemallocator::Jemalloc;

#[macro_use]
extern crate approx;
extern crate arrayvec;
extern crate test;
