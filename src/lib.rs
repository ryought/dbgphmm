#![feature(test)]
#![feature(step_trait)]
#![feature(is_some_and)]
pub mod common;
pub mod distribution;
pub mod e2e;
pub mod float;
pub mod genome;
pub mod graph;
pub mod hashdbg;
pub mod hist;
pub mod hmmv2;
pub mod io;
pub mod kmer;
pub mod min_flow;
pub mod multi_dbg;
pub mod prelude;
pub mod prob;
pub mod random_seq;
pub mod utils;

#[cfg(feature = "python")]
pub mod python;

extern crate jemallocator;
#[global_allocator]
static ALLOC: jemallocator::Jemalloc = jemallocator::Jemalloc;

#[macro_use]
extern crate approx;
extern crate arrayvec;
extern crate test;
