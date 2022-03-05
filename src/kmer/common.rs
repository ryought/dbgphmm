//!
//! Kmer definitions
//!
use crate::common::SeqStyle;

pub trait NullableKmer {
    ///
    /// null <==> NNNNN
    ///
    fn is_null(&self) -> bool;
}

pub trait KmerLike:
    std::marker::Sized + PartialEq + NullableKmer + Eq + std::hash::Hash + Clone + std::fmt::Display
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
        let bases = [b'A', b'C', b'G', b'T'];
        bases
            .iter()
            .map(|&first_base| self.extend_first(first_base))
            .collect()
    }
    ///
    /// XX (k mer) -> [XXA, XXC, XXG, XXT] (k+1 mer)
    ///
    fn succs(&self) -> Vec<Self> {
        let bases = [b'A', b'C', b'G', b'T'];
        bases
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
        self.last() != b'N'
    }
    ///
    /// first base is N
    /// TODO check that the rest is not N
    ///
    fn is_starting(&self) -> bool {
        self.first() == b'N'
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
    // construction
    fn from_bases(bases: &[u8]) -> Self;
    fn to_bases(&self) -> Vec<u8>;
}

///
/// Most fundamental k-mer trait
/// TODO
///
pub trait KmerBase {
    fn k(&self) -> usize;
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
            let mut bases = vec![b'N'; n_prefix];
            bases.extend_from_slice(&self.seq[..n_body]);
            self.index_prefix += 1;
            Some(K::from_bases(&bases))
        } else if self.index <= l - k {
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
                bases.extend_from_slice(&vec![b'N'; n_suffix]);
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
                VecKmer::from_bases(b"NNNA"),
                VecKmer::from_bases(b"NNAT"),
                VecKmer::from_bases(b"NATC"),
                VecKmer::from_bases(b"ATCA"),
                VecKmer::from_bases(b"TCAT"),
                VecKmer::from_bases(b"CATC"),
                VecKmer::from_bases(b"ATCG"),
                VecKmer::from_bases(b"TCGN"),
                VecKmer::from_bases(b"CGNN"),
                VecKmer::from_bases(b"GNNN"),
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
}
