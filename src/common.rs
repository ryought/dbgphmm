//!
//!
//!
pub use petgraph::graph::{EdgeIndex, NodeIndex};
use std::str::FromStr;

/// integer copy number (= occurrence on the genome/sequence)
pub type CopyNum = usize;

/// frequency of the node
pub type Freq = f64;

/// Type of DNA sequence
///
/// If the style of sequence (i.e. circular or linear) matters
/// in your use case, please use `StyledSequence`.
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

///
///
///
#[derive(Clone, Copy, Debug, PartialEq)]
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
    ///
    /// the sequence is fragmented
    /// i.e. tail and head is not connected.
    ///
    pub fn is_fragment(&self) -> bool {
        match self {
            SeqStyle::LinearFragment => true,
            _ => false,
        }
    }
}

impl std::fmt::Display for SeqStyle {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            SeqStyle::Circular => write!(f, "C"),
            SeqStyle::Linear => write!(f, "L"),
            SeqStyle::LinearFragment => write!(f, "F"),
        }
    }
}

impl FromStr for SeqStyle {
    type Err = StyledSequenceParseError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "C" => Ok(SeqStyle::Circular),
            "L" => Ok(SeqStyle::Linear),
            "F" => Ok(SeqStyle::LinearFragment),
            _ => Err(StyledSequenceParseError),
        }
    }
}

///
/// Sequence with style specified.
///
#[derive(Clone, Debug, PartialEq)]
pub struct StyledSequence {
    seq: Sequence,
    style: SeqStyle,
}

impl StyledSequence {
    /// Constructor of Styled Sequence.
    pub fn new(seq: Sequence, style: SeqStyle) -> Self {
        StyledSequence { seq, style }
    }
    pub fn seq(&self) -> &Sequence {
        &self.seq
    }
    pub fn style(&self) -> SeqStyle {
        self.style
    }
}

impl std::fmt::Display for StyledSequence {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}:{}", self.style, sequence_to_string(&self.seq))
    }
}

///
/// Error (unit type) in from_str of StyledSequence and SeqStyle
///
#[derive(Clone, Debug)]
pub struct StyledSequenceParseError;

impl FromStr for StyledSequence {
    type Err = StyledSequenceParseError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let segments: Vec<&str> = s.split(':').collect();
        let style = segments[0].parse::<SeqStyle>()?;
        // TODO sanitize bases
        let seq = segments[1].as_bytes().to_vec();
        Ok(StyledSequence { seq, style })
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn seq_style() {
        let s = SeqStyle::Linear;
        println!("{}", s);
        assert_eq!("L", format!("{}", SeqStyle::Linear));
        assert_eq!("C", format!("{}", SeqStyle::Circular));
        assert_eq!("F", format!("{}", SeqStyle::LinearFragment));
        assert_eq!(SeqStyle::from_str("L").unwrap(), SeqStyle::Linear);
        assert_eq!(SeqStyle::from_str("C").unwrap(), SeqStyle::Circular);
        assert_eq!(SeqStyle::from_str("F").unwrap(), SeqStyle::LinearFragment);
        assert!(SeqStyle::from_str("XX").is_err());
        assert!(SeqStyle::from_str("L ").is_err());
    }
    #[test]
    fn styled_sequence() {
        let s1 = StyledSequence::new(b"ATCGAT".to_vec(), SeqStyle::Circular);
        let e1 = format!("{}", s1);
        assert_eq!(e1, "C:ATCGAT");
        assert_eq!(s1, StyledSequence::from_str(&e1).unwrap());

        let s2 = StyledSequence::new(b"CTCGATCG".to_vec(), SeqStyle::Linear);
        let e2 = "L:CTCGATCG".to_string();
        assert_eq!(s2, StyledSequence::from_str(&e2).unwrap());
    }
}
