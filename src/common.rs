//!
//!
//!
pub use petgraph::graph::{EdgeIndex, NodeIndex};

/// integer copy number (= occurrence on the genome/sequence)
pub type CopyNum = usize;

/// frequency of the node
pub type Freq = f64;

/// Type of DNA sequence
pub type Sequence = Vec<u8>;

/// Convert Sequence(Vec<u8>) into &str
/// useful in displaying
pub fn sequence_to_string(seq: &Sequence) -> &str {
    std::str::from_utf8(seq).unwrap()
}

///
/// short-hand of `NodeIndex::new`
///
pub fn ni(index: usize) -> NodeIndex {
    NodeIndex::new(index)
}

///
/// short-hand of `EdgeIndex::new`
///
pub fn ei(index: usize) -> EdgeIndex {
    EdgeIndex::new(index)
}
