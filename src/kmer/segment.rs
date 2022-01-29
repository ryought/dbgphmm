//! Segment (variable k-mer) definition
use super::common::{KmerBase, KmerLike};

///
/// Segment (variable-length k-mer) struct
///
#[derive(Debug, PartialEq, PartialOrd, Eq, Hash, Clone)]
pub struct Segment {
    bases: Vec<u8>,
    k: usize,
}
