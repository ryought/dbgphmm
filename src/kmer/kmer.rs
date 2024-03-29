//!
//! kmer base struct definitions
//!
pub use super::common::{
    linear_fragment_sequence_to_kmers, linear_sequence_to_kmers as sequence_to_kmers,
    styled_sequence_to_kmers, KmerLike,
};
pub use super::veckmer::VecKmer as Kmer;
use crate::common::NULL_BASE;

use std::collections::HashMap;

pub fn count(seq: &[u8], k: usize) -> HashMap<Kmer, usize> {
    let mut h: HashMap<Kmer, usize> = HashMap::new();
    for i in 0..=&seq.len() - k {
        let kmer = Kmer::from_bases(&seq[i..i + k]);
        // let kmer = &seq[i..i + k];
        // let count = h.entry(kmer).or_insert(0);
        //*count += 1;
        match h.get(&kmer) {
            Some(&occ) => {
                h.insert(kmer, occ + 1);
            }
            _ => {
                h.insert(kmer, 1);
            }
        };
    }
    h
}

/// return X*k
pub fn null_kmer(k: usize) -> Kmer {
    let v = vec![NULL_BASE; k];
    Kmer(v)
}

/// ATTC -> [XXXA, XXAT, XATT]
pub fn starting_kmers(kmer: &Kmer) -> Vec<Kmer> {
    let k = kmer.len();
    let blanks = std::iter::repeat(NULL_BASE).take(k).collect::<Vec<u8>>();
    let heads: Vec<Kmer> = (1..k)
        .map(|i| Kmer([&blanks[i..], &kmer.0[..i]].concat()))
        .collect();
    heads
}

/// ATTC -> [TTCX, TCXX, CXXX]
pub fn ending_kmers(kmer: &Kmer) -> Vec<Kmer> {
    let k = kmer.len();
    let blanks = std::iter::repeat(NULL_BASE).take(k).collect::<Vec<u8>>();
    let tails: Vec<Kmer> = (1..k)
        .map(|i| Kmer([&kmer.0[i..], &blanks[..i]].concat()))
        .collect();
    tails
}

/// ATCGATTC -> Iterator on [NNA, NAT, ATC, TCG, ..., TTC, TCN, CNN]
/// TODO avoid reallocation?
pub fn linear_seq_to_kmers(seq: &[u8], k: usize) -> impl std::iter::Iterator<Item = Kmer> + '_ {
    if seq.len() < k {
        panic!("added sed length={} shorter than k={}", seq.len(), k);
    }
    let first_kmer = Kmer::from_bases(&seq[..k]);
    let last_kmer = Kmer::from_bases(&seq[seq.len() - k..]);
    starting_kmers(&first_kmer)
        .into_iter()
        .chain(seq.windows(k).map(|subseq| Kmer::from_bases(subseq)))
        .chain(ending_kmers(&last_kmer).into_iter())
    /*
    let mut extended_seq = vec![NULL_BASE; k - 1];
    extended_seq.extend_from_slice(seq);
    extended_seq.extend_from_slice(&seq[..k - 1]);
    extended_seq.windows(k).map(|subseq| Kmer::from(subseq))
    */
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kmer_equality() {
        let a = Kmer::from_bases(b"ATCGATTAG");
        let b = Kmer::from_bases(b"ATCGATTAG");
        assert_eq!(a, b);
    }
    #[test]
    fn kmer_adjacency() {
        let a = Kmer::from_bases(b"ATCGATTAG");
        let b = Kmer::from_bases(b"TCGATTAGA");
        assert!(a.adjacent(&b));
    }
    #[test]
    fn kmer_adjacency_fail() {
        let a = Kmer::from_bases(b"ATCGATTAG");
        let b = Kmer::from_bases(b"TCGATTAAA");
        assert!(!a.adjacent(&b));
    }
    #[test]
    fn kmer_last() {
        let a = Kmer::from_bases(b"ATCGATTAG");
        assert_eq!(a.last(), b'G');
    }
    #[test]
    #[should_panic]
    fn k_zero_mer() {
        let a = Kmer::from_bases(b"");
        let x = a.last();
    }
    #[test]
    fn starting() {
        let a = Kmer::from_bases(b"ATCGATTAG");
        assert!(!a.is_starting());
        let a = Kmer::from_bases(b"nnnnAAAAA");
        assert!(!a.is_starting());
        let a = Kmer::from_bases(b"nAGTAAAAA");
        assert!(a.is_starting());
    }
    #[test]
    fn kmer_overlap_pass() {
        {
            let a = Kmer::from_bases(b"ATCG");
            let b = Kmer::from_bases(b"TCGTTA");
            let c = a.overlap(&b, 3);
            println!("{}", c);
            assert_eq!(c, Kmer::from_bases(b"ATCGTTA"));
        }
        {
            let a = Kmer::from_bases(b"ATCG");
            let b = Kmer::from_bases(b"CGTTA");
            let c = a.overlap(&b, 2);
            println!("{}", c);
            assert_eq!(c, Kmer::from_bases(b"ATCGTTA"));
        }
    }
    #[test]
    #[should_panic]
    fn kmer_overlap_fail() {
        {
            let a = Kmer::from_bases(b"ATCG");
            let b = Kmer::from_bases(b"TCGTTA");
            let c = a.overlap(&b, 2);
        }
    }
}
