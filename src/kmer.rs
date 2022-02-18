pub mod common;
pub mod counter;
pub mod kmer;
pub mod quadarray;
pub mod segment;
// pub mod tinykmer;
pub mod veckmer;
// pub mod kmer_with_size;

pub use common::{KmerBase, KmerLike, NullableKmer};
pub use veckmer::VecKmer;
