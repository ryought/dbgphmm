#![feature(test)]

pub mod benchmark;
pub mod cli;
pub mod compressed_dbg;
pub mod cycles;
pub mod dbg;
pub mod distribution;
pub mod graph;
pub mod hmm;
pub mod io;
pub mod kmer;
pub mod optimizer;
pub mod prob;
pub mod random_seq;
pub mod sparse;
pub mod stats;

#[macro_use]
extern crate approx;
extern crate arrayvec;
extern crate test;
