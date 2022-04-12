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

/// Type of Genome, the collection of sequences.
pub type Genome = Vec<Sequence>;

/// Struct for storing multiple emissions, reads.
///
#[derive(Debug, Clone)]
pub struct Reads {
    pub reads: Vec<Sequence>,
}

impl Reads {
    /// Constructor of reads
    pub fn from(reads: Vec<Sequence>) -> Self {
        Reads { reads }
    }
    /// get an iterator over the reads
    pub fn iter(&self) -> impl Iterator<Item = &Sequence> + '_ {
        self.reads.iter()
    }
    /// the number of reads.
    pub fn len(&self) -> usize {
        self.reads.len()
    }
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

///
///
///
#[derive(Clone, Copy, Debug)]
pub enum SeqStyle {
    /// circular sequence
    Circular,
    /// circularized linear sequence
    Linear,
    /// naive linear sequence which is not circular
    LinearFragment,
}

impl SeqStyle {
    /// check if this style is circular or not.
    pub fn is_circular(&self) -> bool {
        match self {
            SeqStyle::Circular => true,
            _ => false,
        }
    }
    pub fn has_prefix(&self) -> bool {
        match self {
            SeqStyle::Linear => true,
            _ => false,
        }
    }
    pub fn has_suffix(&self) -> bool {
        match self {
            SeqStyle::Linear | SeqStyle::Circular => true,
            _ => false,
        }
    }
}
