//!
//! Kmer definitions
//!
use crate::common::{SeqStyle, StyledSequence, NULL_BASE, VALID_BASES};

pub trait NullableKmer {
    ///
    /// Create a null kmer with specified k
    ///
    fn null_kmer(k: usize) -> Self;
    ///
    /// check if this is null-kmer <==> NNNNN
    ///
    fn is_null(&self) -> bool;
    ///
    /// check if the kmer has `N`
    /// e.g. `NNATC, NATTC`.
    ///
    fn has_null(&self) -> bool;
}

pub trait KmerLike:
    std::marker::Sized
    + PartialEq
    + PartialOrd
    + Ord
    + NullableKmer
    + Eq
    + std::hash::Hash
    + Clone
    + std::fmt::Display
    + std::fmt::Debug
{
    ///
    /// k of the k-mer
    ///
    fn len(&self) -> usize;
    ///
    /// k of the k-mer
    /// (an alias of KmerLike.len)
    ///
    fn k(&self) -> usize {
        self.len()
    }
    ///
    /// ABBBB -> A
    ///
    fn first(&self) -> u8;
    ///
    /// AAAAB -> B
    ///
    fn last(&self) -> u8;
    ///
    /// prefix of the kmer
    /// AAAAB -> AAAA
    ///
    fn prefix(&self) -> Self;
    ///
    /// suffix of the kmer
    /// ABBBB -> BBBB
    ///
    fn suffix(&self) -> Self;
    ///
    /// ABBBB and BBBBC is adjacent
    ///
    fn adjacent(&self, other: &Self) -> bool {
        self.suffix() == other.prefix()
    }
    ///
    /// XYYYY -> [YYYYA, YYYYC, YYYYG, YYYYT]
    ///
    fn childs(&self) -> Vec<Self>;
    ///
    /// YYYYZ -> [AYYYY, CYYYY, GYYYY, TYYYY]
    ///
    fn parents(&self) -> Vec<Self>;
    ///
    /// siblings is childs's parents
    /// XYYYY -> [AYYYY, CYYYY, GYYYY, TYYYY]
    ///
    fn siblings(&self) -> Vec<Self>;
    ///
    /// union of childs and parents
    /// XYYYZ -> [
    ///            YYYZA, YYYZC, YYYZG, YYYZT, (childs)
    ///            AXYYY, CXYYY, GXYYY, TXYYY  (parents)
    ///          ]
    ///
    fn neighbors(&self) -> Vec<Self> {
        let mut childs = self.childs();
        let parents = self.parents();
        childs.extend(parents);
        childs
    }
    ///
    /// XX (k mer) -> [AXX, CXX, GXX, TXX] (k+1 mer)
    ///
    fn preds(&self) -> Vec<Self> {
        VALID_BASES
            .iter()
            .map(|&first_base| self.extend_first(first_base))
            .collect()
    }
    ///
    /// XX (k mer) -> [XXA, XXC, XXG, XXT] (k+1 mer)
    ///
    fn succs(&self) -> Vec<Self> {
        VALID_BASES
            .iter()
            .map(|&last_base| self.extend_last(last_base))
            .collect()
    }
    ///
    /// head <==> NNNNX
    ///
    fn is_head(&self) -> bool {
        self.prefix().is_null()
    }
    ///
    /// tail <==> XNNNN
    ///
    fn is_tail(&self) -> bool {
        self.suffix().is_null()
    }
    ///
    /// last base is not N
    ///
    fn is_emitable(&self) -> bool {
        self.last() != NULL_BASE
    }
    ///
    /// first base is N
    /// TODO check that the rest is not N
    ///
    fn is_starting(&self) -> bool {
        self.first() == NULL_BASE
    }
    ///
    /// (YYY, X) -> XYYY
    ///
    fn extend_first(&self, first_base: u8) -> Self;
    ///
    /// (YYY, Z) -> YYYZ
    ///
    fn extend_last(&self, last_base: u8) -> Self;
    ///
    /// (XYYY, YYYZ) (two k mers) -> XYYYZ (k+1 mer)
    ///
    fn join(&self, other: &Self) -> Self {
        if !self.adjacent(other) {
            panic!();
        }
        self.extend_last(other.last())
    }
    ///
    /// upgrade k-mer head into k+1-mer
    ///
    /// NNX -> NNNX
    /// k-mer  k+1-mer
    ///
    fn extend_head(&self) -> Self {
        assert!(self.is_head());
        self.extend_first(NULL_BASE)
    }
    ///
    /// upgrade k-mer tail into k+1-mer
    ///
    /// XNN -> XNNN
    /// k-mer  k+1-mer
    ///
    fn extend_tail(&self) -> Self {
        assert!(self.is_tail());
        self.extend_last(NULL_BASE)
    }
    // construction
    fn from_bases(bases: &[u8]) -> Self;
    fn to_bases(&self) -> Vec<u8>;
}

//
// Sequence <-> Kmers conversion
//

/// Convert linear sequence to a list of kmers
pub fn linear_sequence_to_kmers<'a, K: KmerLike>(
    seq: &'a [u8],
    k: usize,
) -> MarginKmerIterator<'a, K> {
    sequence_to_kmers(seq, k, SeqStyle::Linear)
}

/// Convert linear fragment sequence to a list of kmers
pub fn linear_fragment_sequence_to_kmers<'a, K: KmerLike>(
    seq: &'a [u8],
    k: usize,
) -> MarginKmerIterator<'a, K> {
    sequence_to_kmers(seq, k, SeqStyle::LinearFragment)
}

/// Convert circular sequence to a list of kmers
pub fn circular_sequence_to_kmers<'a, K: KmerLike>(
    seq: &'a [u8],
    k: usize,
) -> MarginKmerIterator<'a, K> {
    sequence_to_kmers(seq, k, SeqStyle::Circular)
}

/// Convert circular or linear sequence to a list of kmers
///
/// # Known Bugs
///
/// * it does not works when `length of seq > k`.
///
fn sequence_to_kmers<'a, K: KmerLike>(
    seq: &'a [u8],
    k: usize,
    seq_style: SeqStyle,
) -> MarginKmerIterator<'a, K> {
    MarginKmerIterator {
        k,
        seq_style,
        index_prefix: 0,
        index_suffix: 0,
        index: 0,
        seq,
        ph: std::marker::PhantomData::<K>,
    }
}

/// Convert circular or linear sequence to a list of kmers
///
/// # Known Bugs
///
/// * it does not works when `length of seq > k`.
///
pub fn styled_sequence_to_kmers<'a, K: KmerLike>(
    s: &'a StyledSequence,
    k: usize,
) -> MarginKmerIterator<'a, K> {
    sequence_to_kmers(s.seq(), k, s.style())
}

///
/// 0 <= index_prefix < k   ('n' * (k-index_prefix)) + seq[:index_prefix]
/// 0 <= index < L - k      seq[index:index+k]
/// 0 <= index_prefix < k   seq[L-index_suffix:L] + ('n' * index_suffix)
///
pub struct MarginKmerIterator<'a, K: KmerLike> {
    k: usize,
    seq_style: SeqStyle,
    index_prefix: usize,
    index_suffix: usize,
    index: usize,
    seq: &'a [u8],
    ph: std::marker::PhantomData<K>,
}

impl<'a, K: KmerLike> Iterator for MarginKmerIterator<'a, K> {
    type Item = K;
    fn next(&mut self) -> Option<K> {
        let k = self.k;
        let l = self.seq.len();
        if self.seq_style.has_prefix() && self.index_prefix < k - 1 {
            // NNNTTT if linear
            let n_prefix = k - 1 - self.index_prefix;
            let n_body = k - n_prefix;
            let mut bases = vec![NULL_BASE; n_prefix];
            bases.extend_from_slice(&self.seq[..n_body]);
            self.index_prefix += 1;
            Some(K::from_bases(&bases))
        } else if l >= k && self.index <= l - k {
            // TTTTTT
            let start = self.index;
            let end = self.index + k;
            let bases = self.seq[start..end].to_vec();
            self.index += 1;
            Some(K::from_bases(&bases))
        } else if self.seq_style.has_suffix() && self.index_suffix < k - 1 {
            // TTTNNN if linear
            // TTTSSS if circular
            let n_suffix = self.index_suffix + 1;
            let n_body = k - n_suffix;
            let mut bases = self.seq[l - n_body..].to_vec();
            if self.seq_style.is_circular() {
                bases.extend_from_slice(&self.seq[..n_suffix]);
            } else {
                bases.extend_from_slice(&vec![NULL_BASE; n_suffix]);
            }
            self.index_suffix += 1;
            Some(K::from_bases(&bases))
        } else {
            // end
            None
        }
    }
}

//
// Tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::veckmer::VecKmer;
    #[test]
    fn seq_to_kmers() {
        let seq = b"ATCATCG";
        println!("linear");
        for kmer in linear_sequence_to_kmers::<VecKmer>(seq, 4) {
            println!("{}", kmer);
        }
        let kmers: Vec<VecKmer> = linear_sequence_to_kmers::<VecKmer>(seq, 4).collect();
        assert_eq!(
            kmers,
            vec![
                VecKmer::from_bases(b"nnnA"),
                VecKmer::from_bases(b"nnAT"),
                VecKmer::from_bases(b"nATC"),
                VecKmer::from_bases(b"ATCA"),
                VecKmer::from_bases(b"TCAT"),
                VecKmer::from_bases(b"CATC"),
                VecKmer::from_bases(b"ATCG"),
                VecKmer::from_bases(b"TCGn"),
                VecKmer::from_bases(b"CGnn"),
                VecKmer::from_bases(b"Gnnn"),
            ]
        );

        let seq = b"ATCATCG";
        println!("circular");
        for kmer in circular_sequence_to_kmers::<VecKmer>(seq, 4) {
            println!("{}", kmer);
        }
        let kmers: Vec<VecKmer> = circular_sequence_to_kmers::<VecKmer>(seq, 4).collect();
        assert_eq!(
            kmers,
            vec![
                VecKmer::from_bases(b"ATCA"),
                VecKmer::from_bases(b"TCAT"),
                VecKmer::from_bases(b"CATC"),
                VecKmer::from_bases(b"ATCG"),
                VecKmer::from_bases(b"TCGA"),
                VecKmer::from_bases(b"CGAT"),
                VecKmer::from_bases(b"GATC"),
            ]
        );

        let seq = b"ATCATCG";
        println!("linear fragment");
        let kmers: Vec<VecKmer> = linear_fragment_sequence_to_kmers::<VecKmer>(seq, 4).collect();
        for kmer in kmers.iter() {
            println!("{}", kmer);
        }
        assert_eq!(
            kmers,
            vec![
                VecKmer::from_bases(b"ATCA"),
                VecKmer::from_bases(b"TCAT"),
                VecKmer::from_bases(b"CATC"),
                VecKmer::from_bases(b"ATCG"),
            ]
        );
    }
    #[test]
    fn seq_to_kmers_edge_cases() {
        // edge case #1: shorter than k
        let seq = b"ATC";
        let k = 4;

        println!("l");
        for kmer in linear_sequence_to_kmers::<VecKmer>(seq, k) {
            println!("{}", kmer);
        }
        let kmers: Vec<VecKmer> = linear_sequence_to_kmers::<VecKmer>(seq, k).collect();
        assert_eq!(
            kmers,
            vec![
                VecKmer::from_bases(b"nnnA"),
                VecKmer::from_bases(b"nnAT"),
                VecKmer::from_bases(b"nATC"),
                VecKmer::from_bases(b"ATCn"),
                VecKmer::from_bases(b"TCnn"),
                VecKmer::from_bases(b"Cnnn"),
            ]
        );

        println!("c");
        for kmer in circular_sequence_to_kmers::<VecKmer>(seq, k) {
            println!("{}", kmer);
        }
        let kmers: Vec<VecKmer> = circular_sequence_to_kmers::<VecKmer>(seq, k).collect();
        assert_eq!(
            kmers,
            vec![
                VecKmer::from_bases(b"ATCA"),
                VecKmer::from_bases(b"TCAT"),
                VecKmer::from_bases(b"CATC"),
            ]
        );

        println!("f");
        for kmer in linear_fragment_sequence_to_kmers::<VecKmer>(seq, k) {
            println!("{}", kmer);
        }
        let kmers: Vec<VecKmer> = linear_fragment_sequence_to_kmers::<VecKmer>(seq, k).collect();
        assert_eq!(kmers.len(), 0);
    }
    #[test]
    fn kmer_extend() {
        let a = VecKmer::from_bases(b"ATCA");
        assert_eq!(a.extend_first(b'A'), VecKmer::from_bases(b"AATCA"));
        assert_eq!(a.extend_last(b'G'), VecKmer::from_bases(b"ATCAG"));

        let a = VecKmer::from_bases(b"nnnA");
        assert!(a.is_head());
        assert_eq!(a.k(), 4);
        let b = a.extend_head();
        assert!(b.is_head());
        assert!(!b.is_tail());
        assert_eq!(b, VecKmer::from_bases(b"nnnnA"));
        assert_eq!(b.k(), 5);

        let a = VecKmer::from_bases(b"Annn");
        assert!(a.is_tail());
        assert_eq!(a.k(), 4);
        let b = a.extend_tail();
        assert!(b.is_tail());
        assert_eq!(b, VecKmer::from_bases(b"Annnn"));
        assert_eq!(b.k(), 5);
    }
    #[test]
    fn kmer_null() {
        let a = VecKmer::from_bases(b"nnnn");
        assert!(a.is_null());
        assert!(a.has_null());

        let b = VecKmer::from_bases(b"nnnT");
        assert!(!b.is_null());
        assert!(b.has_null());

        let c = VecKmer::from_bases(b"TGAC");
        assert!(!c.is_null());
        assert!(!c.has_null());
    }
    #[test]
    fn kmer_order() {
        let mut kmers = vec![
            VecKmer::from_bases(b"ATCG"),
            VecKmer::from_bases(b"AAAA"),
            VecKmer::from_bases(b"nnnn"),
            VecKmer::from_bases(b"TTTT"),
            VecKmer::from_bases(b"CTAG"),
        ];
        kmers.sort();
        for kmer in kmers.iter() {
            println!("{}", kmer);
        }
    }
}
