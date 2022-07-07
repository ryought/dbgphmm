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
//! ## Types of data generations
//! defined in em::e2e::fragments and full_length
//!
//! * FullLengthRead: `full_length`
//! * FragmentRead: `fragments`
//!
pub mod compression;
mod fragments;
mod full_length;
mod runner;
mod tandem_repeat;
mod tandem_repeat_v3;
