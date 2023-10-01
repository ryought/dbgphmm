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
use crate::common::CopyNum;
pub use crate::genome::Genome;
use crate::graph::genome_graph::{GenomeGraphPos, GenomeGraphPosVec};
use crate::hmmv2::hint::Mapping;
use crate::kmer::kmer::KmerLike;
use fnv::FnvHashMap as HashMap;
use itertools::Itertools;
use petgraph::graph::NodeIndex;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_with::{
    serde_as, DeserializeAs, DeserializeFromStr, DisplayFromStr, SerializeAs, SerializeDisplay,
};
use std::io::Write;
use std::str::FromStr;

#[cfg(feature = "python")]
use pyo3::prelude::*;

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
pub trait Seq: AsRef<Bases> + Sync + Clone {
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
impl<T: AsRef<Bases> + Sync + Clone> Seq for T {}
// TODO
// impl<T: Seq> std::fmt::Display for T {
//     fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
//         write!(f, "hoge")
//     }
// }

/// Convert Sequence(Vec<u8>) into &str
/// useful in displaying
pub fn sequence_to_string<T: AsRef<Bases>>(seq: &T) -> &str {
    std::str::from_utf8(seq.as_ref()).unwrap()
}

///
///
///
pub fn starts_and_ends_of_genome<K: KmerLike>(
    genome: &Genome,
    k: usize,
) -> (Vec<(K, CopyNum)>, Vec<(K, CopyNum)>) {
    // mapping from kmer to copynum
    let mut starts: HashMap<K, CopyNum> = HashMap::default();
    let mut ends: HashMap<K, CopyNum> = HashMap::default();
    for hap in genome {
        let (starts_hap, ends_hap) = hap.start_and_end_kmer::<K>(k);
        for start in starts_hap.into_iter() {
            *starts.entry(start).or_insert(0) += 1;
        }
        for end in ends_hap.into_iter() {
            *ends.entry(end).or_insert(0) += 1;
        }
    }
    // convert HashMap<K, CopyNum> into Vec<(K, CopyNum)> (key value list)
    (starts.into_iter().collect(), ends.into_iter().collect())
}

///
/// Struct for storing multiple emissions, reads.
///
/// Read type S in ReadCollection<S> can be
/// * Vec<u8>
/// * PositionedSequence
/// * (TODO StyledSequence)
///
#[serde_as]
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(bound = "S: StoreableTypeTrait")]
pub struct ReadCollection<S: Seq> {
    #[serde_as(as = "Vec<StoreableType>")]
    pub reads: Vec<S>,
}

//
// support serialize/deserialize of Vec<u8> as Bases
//
trait StoreableTypeTrait {
    fn to_store_string(&self) -> String;
    fn from_store_str(s: &str) -> Self;
}
impl StoreableTypeTrait for Vec<u8> {
    fn to_store_string(&self) -> String {
        self.to_str().to_owned()
    }
    fn from_store_str(s: &str) -> Vec<u8> {
        s.to_string().into_bytes()
    }
}
pub struct StoreableType;
impl<T: StoreableTypeTrait> SerializeAs<T> for StoreableType {
    fn serialize_as<S>(source: &T, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(&source.to_store_string())
    }
}
impl<'de, T: StoreableTypeTrait> DeserializeAs<'de, T> for StoreableType {
    fn deserialize_as<D>(deserializer: D) -> Result<T, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = String::deserialize(deserializer).map_err(serde::de::Error::custom)?;
        Ok(T::from_store_str(&s))
    }
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
    /// total number of bases of all reads
    pub fn total_bases(&self) -> usize {
        self.reads.iter().map(|read| read.as_ref().len()).sum()
    }
    /// average length
    pub fn average_length(&self) -> usize {
        self.total_bases() / self.len()
    }
    pub fn coverage(&self, genome_size: usize) -> f64 {
        self.total_bases() as f64 / genome_size as f64
    }
    /// show reads
    ///
    /// ```text
    /// read#1   ATCGTAGCT
    /// read#2   ATCGTA
    /// ```
    ///
    pub fn show_reads(&self) {
        for (i, read) in self.iter().enumerate() {
            println!("# read#{}\t{}", i, read.to_str());
        }
    }
    ///
    /// Dump FASTA file of reads
    /// Position data is ignored
    ///
    pub fn to_fasta<P: AsRef<std::path::Path>>(&self, path: P) -> std::io::Result<()> {
        let file = std::fs::File::create(path).unwrap();
        let mut writer = bio::io::fasta::Writer::new(file);

        for (i, r) in self.iter().enumerate() {
            writer.write(&format!("r{}", i), None, r.as_ref())?;
        }

        Ok(())
    }
}

pub fn sanitize_bases(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .enumerate()
        .map(|(i, base)| match base {
            b'A' | b'a' => b'A',
            b'C' | b'c' => b'C',
            b'G' | b'g' => b'G',
            b'T' | b't' => b'T',
            &c => {
                panic!("Non DNA base '{}' appeared in seq", c as char)
            }
        })
        .collect()
}

impl ReadCollection<Sequence> {
    ///
    /// Load FASTA file
    ///
    pub fn from_fasta<P: AsRef<std::path::Path>>(path: P) -> std::io::Result<Self> {
        let file = std::fs::File::open(path).unwrap();
        let reader = bio::io::fasta::Reader::new(file);

        let reads: Vec<Vec<u8>> = reader
            .records()
            .map(|r| sanitize_bases(r.unwrap().seq()))
            .collect();

        Ok(ReadCollection { reads })
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
    pub fn to_reads(self) -> Reads {
        let reads: Vec<Sequence> = self.reads.into_iter().map(|pos_seq| pos_seq.seq).collect();
        Reads::from(reads)
    }
    ///
    /// make reads' strands same
    ///
    pub fn justify_strand(self) -> PositionedReads {
        let reads: Vec<_> = self
            .reads
            .into_iter()
            .map(|pos_seq| {
                if pos_seq.is_revcomp() {
                    pos_seq.revcomp()
                } else {
                    pos_seq
                }
            })
            .collect();
        PositionedReads::from(reads)
    }
    ///
    /// Dump SAM file of reads
    ///
    pub fn to_sam<P: AsRef<std::path::Path>>(&self, path: P) -> std::io::Result<()> {
        let mut file = std::fs::File::create(path).unwrap();

        for (i, r) in self.iter().enumerate() {
            writeln!(file, "{}", r.to_sam_string(&format!("r{}", i)))?;
        }

        Ok(())
    }
}

impl<'a, S: Seq> IntoIterator for &'a ReadCollection<S> {
    type Item = &'a S;
    type IntoIter = std::slice::Iter<'a, S>;
    fn into_iter(self) -> std::slice::Iter<'a, S> {
        self.reads.iter()
    }
}

impl<'a, S: Seq> IntoParallelIterator for &'a ReadCollection<S> {
    type Item = &'a S;
    type Iter = rayon::slice::Iter<'a, S>;
    fn into_par_iter(self) -> rayon::slice::Iter<'a, S> {
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
#[derive(Clone, Debug, PartialEq, SerializeDisplay, DeserializeFromStr)]
pub struct StyledSequence {
    pub seq: Sequence,
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
    pub fn start_and_end_kmer<K: KmerLike>(&self, k: usize) -> (Vec<K>, Vec<K>) {
        match self.style() {
            SeqStyle::Linear => {
                let n = self.len();
                (
                    vec![K::from_bases(&self.seq[0..k])],
                    vec![K::from_bases(&self.seq[n - k..n])],
                )
            }
            _ => (Vec::new(), Vec::new()),
        }
    }
}

// impl StyledSequence {
//     #[new]
//     fn __new__(s: &str) -> Self {
//         StyledSequence::from_str(s).unwrap()
//     }
//     fn __repr__(&self) -> String {
//         self.to_string()
//     }
// }

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
// required for serde
impl std::fmt::Display for StyledSequenceParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "styleparseeerror")
    }
}

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
/// supported traits
/// * Display, FromStr
/// * ScoreableTypeTrait (using Display and FromStr)
/// * serde::Serialize, Deserialize (using ScoreableTypeTrait)
/// * AsRef<Bases>, Seq (by ignoring position information)
///
#[derive(Clone, Debug, PartialEq, SerializeDisplay, DeserializeFromStr)]
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
        assert_eq!(seq.len(), origins.len());
        PositionedSequence {
            seq,
            origins,
            is_revcomp,
        }
    }
    pub fn len(&self) -> usize {
        self.seq.len()
    }
    pub fn origins(&self) -> &GenomeGraphPosVec {
        &self.origins
    }
    pub fn seq(&self) -> &Sequence {
        &self.seq
    }
    /// origin position of first base
    pub fn head_origin(&self) -> GenomeGraphPos {
        *self.origins.first().unwrap()
    }
    pub fn first_match_origin(&self) -> GenomeGraphPos {
        *self
            .origins
            .iter()
            .find(|origin| origin.is_match())
            .expect("sequence has no valid origin info")
    }
    ///
    /// origin node (NodeIndex) of this positioned sequence.
    /// Determined by the first (Match) origin of the read.
    ///
    pub fn origin_node(&self) -> NodeIndex {
        self.first_match_origin().node().unwrap()
    }
    ///
    /// origin position (location of bases; usize) of this positioned sequence.
    /// Determined by the first (Match) origin of the read.
    ///
    pub fn origin_pos(&self) -> usize {
        self.first_match_origin().pos().unwrap()
    }
    /// origin position of last base
    pub fn tail_origin(&self) -> GenomeGraphPos {
        *self.origins.last().unwrap()
    }
    pub fn is_revcomp(&self) -> bool {
        self.is_revcomp
    }
    /// Make seq reverse complement
    pub fn revcomp(self) -> Self {
        let mut origins = self.origins;
        origins.reverse();
        PositionedSequence {
            seq: self.seq.to_revcomp(),
            origins,
            is_revcomp: !self.is_revcomp,
        }
    }
}

#[derive(Clone, Debug, Copy, PartialEq)]
enum Aln {
    M,
    I,
    D,
}

/// Run Length Encoding
/// [A,A,A,B,B,C,B,B] -> [(A,3),(B,2),(C,1),(B,2)]
fn run_length_encoding<T: PartialEq + Copy>(xs: &[T]) -> Vec<(T, usize)> {
    let mut ret = Vec::new();

    // if xs.len() == 0
    if xs.is_empty() {
        return ret;
    }

    // assuming xs.len() >= 1
    let mut x0 = xs[0];
    let mut n = 1;
    for &x in &xs[1..] {
        if x == x0 {
            // continue
            n += 1;
        } else {
            // different symbol
            ret.push((x0, n));
            x0 = x;
            n = 1;
        }
    }
    ret.push((x0, n));
    ret
}

///
/// Output-related methods of PositionedSequence
///
impl PositionedSequence {
    /// PAF format
    ///
    /// https://lh3.github.io/minimap2/minimap2.html#10
    pub fn to_paf_string(&self) -> String {
        unimplemented!();
    }
    /// create `Vec<Aln> = [M/I/D]`
    fn alignments(&self) -> Vec<Aln> {
        let mut pos = None;
        let mut ret = Vec::new();

        for &origin in self.origins.iter() {
            if origin.is_match() {
                let x = origin.pos().unwrap();
                if let Some(y) = pos {
                    let n_del = x - y - 1;
                    for _ in 0..n_del {
                        ret.push(Aln::D);
                    }
                }
                ret.push(Aln::M);
                pos = Some(x);
            } else {
                ret.push(Aln::I);
            }
        }
        ret
    }
    /// generate cigar
    /// that is run-length encoding of `self.alignments()`
    fn cigar(&self) -> String {
        assert!(!self.is_revcomp);
        run_length_encoding(&self.alignments())
            .into_iter()
            .map(|(aln, count)| match aln {
                Aln::M => format!("M{}", count),
                Aln::I => format!("I{}", count),
                Aln::D => format!("D{}", count),
            })
            .join("")
    }
    /// SAM format
    ///
    /// https://samtools.github.io/hts-specs/SAMv1.pdf
    fn to_sam_string(&self, name: &str) -> String {
        let node = self.origin_node();
        let pos = self.origin_pos();
        format!(
            "{}\t0\tg{}\t{}\t255\t{}\t*\t0\t0\t{}\t*",
            name,
            node.index(),
            pos + 1,
            self.cigar(),
            self.seq().to_str(),
        )
    }
    ///
    /// To aligned two-row string representation
    ///
    /// ```text
    /// ATGCGA-CGTGG
    ///   G     C
    /// ```
    ///
    /// Deletion: `-` in 1st-row
    /// Insertion: `X` in 2nd-row
    ///
    pub fn to_aligned_str(&self) -> [String; 2] {
        // // let n = self.origins.iter().filter(|origin| origin.is_match())
        // let mut row1 = String::new();
        // let mut row2 = String::new();
        // for i in 0..self.len() {
        //     if self.origins[i].is_match() {
        //         row1.push(self.seq[i] as char);
        //         row2.push(' ');
        //     } else {
        //         // row1.push(self.seq[i] as char);
        //         row2.push(self.seq[i] as char);
        //     }
        // }
        // [row1, row2]
        unimplemented!();
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
            "{}:{}:{}",
            self.seq.to_str(),
            if self.is_revcomp() { '-' } else { '+' },
            self.origins
                .iter()
                .map(|origin| origin.to_string())
                .join(","),
        )
    }
}
#[derive(Clone, Debug)]
pub struct PositionedSequenceParseError;
impl std::fmt::Display for PositionedSequenceParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "PositionedSequenceParseError")
    }
}
impl FromStr for PositionedSequence {
    type Err = PositionedSequenceParseError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let segments: Vec<&str> = s.split(':').collect();
        if segments.len() == 3 {
            let seq = segments[0].as_bytes().to_vec();
            let is_revcomp = if segments[1] == "+" {
                false
            } else if segments[1] == "-" {
                true
            } else {
                return Err(PositionedSequenceParseError);
            };
            let origins: Vec<GenomeGraphPos> =
                segments[2].split(',').map(|o| o.parse().unwrap()).collect();
            Ok(PositionedSequence {
                seq,
                origins,
                is_revcomp,
            })
        } else {
            Err(PositionedSequenceParseError)
        }
    }
}

impl StoreableTypeTrait for PositionedSequence {
    fn to_store_string(&self) -> String {
        self.to_string()
    }
    fn from_store_str(s: &str) -> Self {
        s.parse().unwrap()
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ni;
    use crate::kmer::common::kmers_to_string;
    use crate::kmer::VecKmer;
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
    #[test]
    fn reads_serialize() {
        {
            // with no hints
            let reads =
                ReadCollection::from(vec![b"ATCGATTCGTA".to_vec(), b"TTTTTTGTGGGGTG".to_vec()]);
            let json = serde_json::to_string(&reads).unwrap();
            println!("{}", json);
            assert_eq!(json, "{\"reads\":[\"ATCGATTCGTA\",\"TTTTTTGTGGGGTG\"]}");
            let reads2: Reads = serde_json::from_str(&json).unwrap();
            println!("{:?}", reads2);
            assert_eq!(reads, reads2);
        }
    }
    #[test]
    fn styled_seq_serialize() {
        let ts = vec![1, 2, 3];
        println!("{}", serde_json::to_string(&ts).unwrap());

        // single styledsequence
        let x0 = StyledSequence::new(b"ATCGAT".to_vec(), SeqStyle::Circular);
        let json = serde_json::to_string(&x0).unwrap();
        assert_eq!(json, "\"C:ATCGAT\"");
        let x: StyledSequence = serde_json::from_str(&json).unwrap();
        println!("x={}", x);
        assert_eq!(x, x0);

        // multiple styledsequences
        let xs0 = vec![
            StyledSequence::new(b"ATCGAT".to_vec(), SeqStyle::Circular),
            StyledSequence::new(b"GGGC".to_vec(), SeqStyle::Linear),
        ];
        let json = serde_json::to_string(&xs0).unwrap();
        println!("{}", json);
        assert_eq!(json, "[\"C:ATCGAT\",\"L:GGGC\"]");
        let xs: Vec<StyledSequence> = serde_json::from_str(&json).unwrap();
        assert_eq!(xs, xs0);
    }
    #[test]
    fn positioned_seq_serialize() {
        let s1 = PositionedSequence::new(
            b"ATCGT".to_vec(),
            vec![
                GenomeGraphPos::new_match(ni(0), 0),
                GenomeGraphPos::new_match(ni(0), 1),
                GenomeGraphPos::new_match(ni(0), 2),
                GenomeGraphPos::new_ins(),
                GenomeGraphPos::new_match(ni(0), 3),
            ],
            true,
        );
        println!("{}", s1);
        let s1b = PositionedSequence::from_str(&s1.to_string()).unwrap();
        assert_eq!(s1, s1b);
    }
    #[test]
    fn starts_and_ends() {
        let genome = Genome::new(vec![
            // hap0
            StyledSequence::new(b"ATCGATTTAGC".to_vec(), SeqStyle::Linear),
            // hap1
            StyledSequence::new(b"GGGCGGCTGCTG".to_vec(), SeqStyle::Linear),
            // hap2
            // start kmer of hap1 and hap2 are same (GGGC).
            StyledSequence::new(b"GGGCGGCTGCTGC".to_vec(), SeqStyle::Linear),
            // hap3
            // circular does not have start/end
            StyledSequence::new(b"GGGCGGCTGCTG".to_vec(), SeqStyle::Circular),
        ]);
        let (starts, ends) = starts_and_ends_of_genome::<VecKmer>(&genome, 4);
        for (kmer, copy_num) in starts.iter() {
            println!("start {} x{}", kmer, copy_num);
        }
        for (kmer, copy_num) in ends.iter() {
            println!("end {} x{}", kmer, copy_num);
        }
        assert_eq!(
            starts,
            vec![
                (VecKmer::from_bases(b"GGGC"), 2),
                (VecKmer::from_bases(b"ATCG"), 1),
            ],
        );
        assert_eq!(
            ends,
            vec![
                (VecKmer::from_bases(b"TAGC"), 1),
                (VecKmer::from_bases(b"CTGC"), 1),
                (VecKmer::from_bases(b"GCTG"), 1),
            ],
        );
    }
    #[test]
    fn positioned_seq_alignment() {
        let s1 = PositionedSequence::new(
            b"ATCGTTCG".to_vec(),
            vec![
                GenomeGraphPos::new_match(ni(0), 0),
                GenomeGraphPos::new_match(ni(0), 1),
                GenomeGraphPos::new_match(ni(0), 3),
                GenomeGraphPos::new_match(ni(0), 4),
                GenomeGraphPos::new_match(ni(0), 5),
                GenomeGraphPos::new_match(ni(0), 8),
                GenomeGraphPos::new_ins(),
                GenomeGraphPos::new_match(ni(0), 9),
            ],
            false,
        );
        println!("{:?}", s1.alignments());
        assert_eq!(
            s1.alignments(),
            vec![
                Aln::M, // x[0]
                Aln::M, // x[1]
                Aln::D, // x[2]
                Aln::M, // x[3]
                Aln::M, // x[4]
                Aln::M, // x[5]
                Aln::D, // x[6]
                Aln::D, // x[7]
                Aln::M, // x[8]
                Aln::I,
                Aln::M, // x[9]
            ]
        );
        assert_eq!(s1.cigar(), "M2D1M3D2M1I1M1");
        println!("{}", s1.to_sam_string("r1"));
        assert_eq!(
            s1.to_sam_string("r1"),
            "r1\t0\tg0\t1\t255\tM2D1M3D2M1I1M1\t*\t0\t0\tATCGTTCG\t*"
        );

        let s2 = PositionedSequence::new(
            b"GTCGTTCGG".to_vec(),
            vec![
                GenomeGraphPos::new_ins(),
                GenomeGraphPos::new_ins(),
                GenomeGraphPos::new_ins(),
                GenomeGraphPos::new_match(ni(1), 10),
                GenomeGraphPos::new_ins(),
                GenomeGraphPos::new_match(ni(1), 15),
                GenomeGraphPos::new_match(ni(1), 17),
                GenomeGraphPos::new_ins(),
                GenomeGraphPos::new_match(ni(1), 18),
            ],
            false,
        );
        assert_eq!(
            s2.alignments(),
            vec![
                Aln::I,
                Aln::I,
                Aln::I,
                Aln::M, // x[10]
                Aln::I,
                Aln::D, // x[11]
                Aln::D, // x[12]
                Aln::D, // x[13]
                Aln::D, // x[14]
                Aln::M, // x[15]
                Aln::D, // x[16]
                Aln::M, // x[17]
                Aln::I,
                Aln::M, // x[18]
            ]
        );
        println!("{:?}", s2.alignments());
        assert_eq!(s2.cigar(), "I3M1I1D4M1D1M1I1M1");
        println!("{}", s2.to_sam_string("r2"));
        assert_eq!(
            s2.to_sam_string("r2"),
            "r2\t0\tg1\t11\t255\tI3M1I1D4M1D1M1I1M1\t*\t0\t0\tGTCGTTCGG\t*"
        );

        // let aligned = s1.to_aligned_str();
        // println!("{}\n{}", aligned[0], aligned[1]);
    }
    #[test]
    fn read_collection_fasta() {
        let r = ReadCollection {
            reads: vec![b"ATCGGATC".to_vec(), b"TTTTTTA".to_vec(), b"GCTAG".to_vec()],
        };
        r.to_fasta("tmp.fa").unwrap();
        r.show_reads();

        let r2 = ReadCollection::from_fasta("tmp.fa").unwrap();
        r2.show_reads();

        assert_eq!(r, r2);
    }
}
