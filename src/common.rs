//!
//!
//!
pub use petgraph::graph::{EdgeIndex, NodeIndex};
pub mod collection;
pub use collection::{
    sequence_to_string, Bases, Genome, PositionedReads, PositionedSequence, Reads, Seq, SeqStyle,
    Sequence, StyledSequence, StyledSequenceParseError,
};

/// integer copy number (= occurrence on the genome/sequence)
pub type CopyNum = usize;

/// frequency of the node
pub type Freq = f64;

///
/// null base
///
pub const NULL_BASE: u8 = b'n';

///
/// Array of valid DNA bases
///
pub const VALID_BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

///
/// Array of all bases
///
pub const BASES: [u8; 5] = [b'A', b'C', b'G', b'T', NULL_BASE];

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
