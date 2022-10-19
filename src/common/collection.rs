//!
//! Sequence and its collections
//!
//! ## Single Sequence
//!
//! * `Sequence`
//! * `StyledSequence`
//! * `PositionedSequence`
//!
//! ## Sequences
//!
//! * `Genome`: simple vector
//! * `Reads`: read collections
//!
use crate::graph::genome_graph::{GenomeGraphPos, GenomeGraphPosVec};
use itertools::Itertools;
use pyo3::prelude::*;
use rayon::prelude::*;
use std::str::FromStr;

//
// Single Sequence
//

/// Type of DNA sequence
///
/// If the style of sequence (i.e. circular or linear) matters
/// in your use case, please use `StyledSequence`.
pub type Sequence = Vec<u8>;

/// Type of Bases as array
///
/// It is used in `AsRef<Bases>` or `&Bases`
pub type Bases = [u8];

///
/// Get the complemented base A <=> T, G <=> C
///
pub fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        n => n,
    }
}

/// Seq trait
/// It can be converted into &Bases with `as_ref()`.
///
/// * `as_ref`
/// * `to_str`
/// * `to_revcomp`
///
pub trait Seq: AsRef<Bases> {
    ///
    /// convert bases into &str for displaying
    ///
    fn to_str(&self) -> &str {
        std::str::from_utf8(self.as_ref()).unwrap()
    }
    ///
    /// convert into reverse complemented sequence
    ///
    fn to_revcomp(&self) -> Sequence {
        self.as_ref()
            .iter()
            .rev()
            .map(|&base| complement(base))
            .collect()
    }
}
impl<T: AsRef<Bases>> Seq for T {}

/// Convert Sequence(Vec<u8>) into &str
/// useful in displaying
pub fn sequence_to_string<T: AsRef<Bases>>(seq: &T) -> &str {
    std::str::from_utf8(seq.as_ref()).unwrap()
}

/// Type of Genome, the collection of sequences.
pub type Genome = Vec<StyledSequence>;

///
/// calculate genome size of the given genome.
///
pub fn genome_size(genome: &Genome) -> usize {
    genome.iter().map(|seq| seq.len()).sum()
}

/// Struct for storing multiple emissions, reads.
///
#[derive(Debug, Clone)]
pub struct ReadCollection<S: Seq> {
    pub reads: Vec<S>,
}

///
/// Reads
///
/// ReadCollection of Sequences
///
pub type Reads = ReadCollection<Sequence>;

///
/// PositionedReads
///
/// ReadCollection of PositionedSequences
///
pub type PositionedReads = ReadCollection<PositionedSequence>;

impl<S: Seq> ReadCollection<S> {
    /// Constructor of reads
    pub fn from(reads: Vec<S>) -> Self {
        ReadCollection { reads }
    }
    /// get an iterator over the reads
    pub fn iter(&self) -> impl Iterator<Item = &S> + '_ {
        self.reads.iter()
    }
    /// the number of reads.
    pub fn len(&self) -> usize {
        self.reads.len()
    }
}

impl<S: Seq> std::ops::Index<usize> for ReadCollection<S> {
    type Output = S;
    fn index(&self, index: usize) -> &Self::Output {
        &self.reads[index]
    }
}

impl PositionedReads {
    ///
    /// Remove position information and convert to `Reads`.
    ///
    /// * `justify_strand` (bool)
    ///     if true, align all reads in forward strand
    ///     by revcomping backward reads.
    ///
    pub fn to_reads(self, justify_strand: bool) -> Reads {
        let reads: Vec<Sequence> = self
            .reads
            .into_iter()
            .map(|pos_seq| {
                if justify_strand && pos_seq.is_revcomp() {
                    pos_seq.seq.to_revcomp()
                } else {
                    pos_seq.seq
                }
            })
            .collect();
        Reads::from(reads)
    }
}

impl<'a> IntoIterator for &'a Reads {
    type Item = &'a Sequence;
    type IntoIter = std::slice::Iter<'a, Sequence>;
    fn into_iter(self) -> std::slice::Iter<'a, Sequence> {
        self.reads.iter()
    }
}

impl<'a> IntoParallelIterator for &'a Reads {
    type Item = &'a Sequence;
    type Iter = rayon::slice::Iter<'a, Sequence>;
    fn into_par_iter(self) -> rayon::slice::Iter<'a, Sequence> {
        self.reads.par_iter()
    }
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
#[pyclass]
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
    /// Construct StyledSequence with SeqStyle::Circular
    pub fn circular(seq: Sequence) -> Self {
        StyledSequence::new(seq, SeqStyle::Circular)
    }
    /// Construct StyledSequence with SeqStyle::Linear
    pub fn linear(seq: Sequence) -> Self {
        StyledSequence::new(seq, SeqStyle::Linear)
    }
    /// Construct StyledSequence with SeqStyle::LinearFragment
    pub fn linear_fragment(seq: Sequence) -> Self {
        StyledSequence::new(seq, SeqStyle::LinearFragment)
    }
    pub fn seq(&self) -> &Sequence {
        &self.seq
    }
    pub fn to_seq(self) -> Sequence {
        self.seq
    }
    pub fn style(&self) -> SeqStyle {
        self.style
    }
    /// length of the sequence
    pub fn len(&self) -> usize {
        self.seq.len()
    }
}

#[pymethods]
impl StyledSequence {
    #[new]
    fn __new__(s: &str) -> Self {
        StyledSequence::from_str(s).unwrap()
    }
    fn __repr__(&self) -> String {
        self.to_string()
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

impl AsRef<Bases> for StyledSequence {
    fn as_ref(&self) -> &Bases {
        &self.seq
    }
}

impl AsRef<StyledSequence> for StyledSequence {
    fn as_ref(&self) -> &StyledSequence {
        self
    }
}

///
/// PositionedSequence
///
#[derive(Clone, Debug)]
pub struct PositionedSequence {
    /// sequence (bases) of the read
    seq: Sequence,
    /// origin position of each base of the read
    origins: GenomeGraphPosVec,
    /// is reverse complemented or not
    is_revcomp: bool,
}

impl PositionedSequence {
    /// Constructor of positioned sequence.
    pub fn new(seq: Sequence, origins: GenomeGraphPosVec, is_revcomp: bool) -> Self {
        PositionedSequence {
            seq,
            origins,
            is_revcomp,
        }
    }
    pub fn origins(&self) -> &GenomeGraphPosVec {
        &self.origins
    }
    /// origin position of first base
    pub fn head_origin(&self) -> GenomeGraphPos {
        *self.origins.first().unwrap()
    }
    /// origin position of last base
    pub fn tail_origin(&self) -> GenomeGraphPos {
        *self.origins.last().unwrap()
    }
    pub fn is_revcomp(&self) -> bool {
        self.is_revcomp
    }
}

impl AsRef<Bases> for PositionedSequence {
    fn as_ref(&self) -> &Bases {
        &self.seq
    }
}

impl std::fmt::Display for PositionedSequence {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{} (revcomp={}, origins={})",
            self.seq.to_str(),
            self.is_revcomp(),
            self.origins
                .iter()
                .map(|origin| origin.to_string())
                .join(","),
        )
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
    #[test]
    fn seq_revcomp() {
        let s1 = b"ATCGGCCC".to_vec();
        let s2 = s1.to_revcomp();
        println!("{}", s1.to_str());
        println!("{}", s2.to_str());
        assert_eq!(s2, b"GGGCCGAT");
    }
}
