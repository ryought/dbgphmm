#![feature(test)]
#![feature(vec_retain_mut)]
#![feature(int_abs_diff)]
#![feature(step_trait)]
pub mod cli;
pub mod common;
pub mod compressed_dbg;
pub mod cycles;
pub mod dbg;
pub mod distribution;
pub mod e2e;
pub mod genome;
pub mod graph;
pub mod greedy;
pub mod hist;
pub mod hmm;
pub mod hmmv2;
pub mod inspect;
pub mod io;
pub mod json;
pub mod kmer;
pub mod mcmc;
pub mod optimizer;
pub mod min_flow;
pub mod playground;
pub mod prelude;
pub mod prob;
pub mod python;
pub mod random_seq;
pub mod reads;
pub mod segment;
pub mod stats;
pub mod utils;
pub mod veclike;
pub mod vector;

extern crate jemallocator;
#[global_allocator]
static ALLOC: jemallocator::Jemalloc = jemallocator::Jemalloc;

#[macro_use]
extern crate approx;
extern crate arrayvec;
extern crate test;
