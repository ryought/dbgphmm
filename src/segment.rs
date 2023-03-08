//!
//! Kmer
//!
use crate::common::Sequence;

///
pub struct Segment {
    seq: Sequence,
}

impl Segment {
    pub fn new(bases: &[u8]) -> Self {
        Segment {
            seq: bases.to_owned(),
        }
    }
    ///
    /// AAAXXX and XXXB is adjacent when k=3
    ///
    pub fn adjacent(&self, other: &Self, k: usize) -> bool {
        true
    }
}

impl std::fmt::Display for Segment {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        unimplemented!();
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn f() {}
}
