//!
//! End-to-end genome inference tests
//!
//! ## Overview
//!
//! It starts from a raw dbg constructed from reads.
//!
//! 1. generate genome
//! 2. generate reads from linear graph
//!
//!
//! ## Types of genomes
//! defined in em::e2e::genome
//! *
//!
//! ## Types of data generations
//! defined in em::e2e::fragments and full_length
//!
//! * FullLengthRead: `full_length`
//! * FragmentRead: `fragments`
//!
mod fragments;
mod full_length;
mod genome;
mod runner;
mod tandem_repeat;
