//! VecKmer definitions
use super::common::{KmerLike, NullableKmer};
use crate::common::{BASES, NULL_BASE, VALID_BASES};
use pyo3::prelude::*;

///
/// Kmer for any k, by using Vec<u8> as a internal struct
///
#[pyclass]
#[derive(Debug, PartialEq, PartialOrd, Ord, Eq, Hash, Clone)]
pub struct VecKmer(pub Vec<u8>);

///
/// short-hand of `VecKmer::from_bases`.
///
pub fn kmer(bases: &[u8]) -> VecKmer {
    VecKmer::from_bases(bases)
}

impl NullableKmer for VecKmer {
    fn null_kmer(k: usize) -> Self {
        VecKmer(vec![NULL_BASE; k])
    }
    fn is_null(&self) -> bool {
        self.0.iter().all(|&x| x == NULL_BASE)
    }
    fn has_null(&self) -> bool {
        self.0.iter().any(|&x| x == NULL_BASE)
    }
}

impl KmerLike for VecKmer {
    fn len(&self) -> usize {
        self.0.len()
    }
    fn adjacent(&self, other: &VecKmer) -> bool {
        let (_, a_suffix) = self.0.split_first().expect("k should be >1");
        let (_, b_prefix) = other.0.split_last().expect("k should be >1");
        a_suffix == b_prefix
    }
    fn first(&self) -> u8 {
        let (first, _) = self.0.split_first().expect("k should be >=1");
        *first
    }
    fn last(&self) -> u8 {
        let (last, _) = self.0.split_last().expect("k should be >=1");
        *last
    }
    fn prefix(&self) -> VecKmer {
        let (_, prefix) = self.0.split_last().expect("k should be >1");
        VecKmer(prefix.to_vec())
    }
    fn suffix(&self) -> VecKmer {
        let (_, suffix) = self.0.split_first().expect("k should be >1");
        VecKmer(suffix.to_vec())
    }
    fn childs(&self) -> Vec<VecKmer> {
        let (_, suffix) = self.0.split_first().unwrap();
        let childs = BASES
            .iter()
            .map(|last_base| {
                let mut v = suffix.to_vec();
                v.push(*last_base);
                VecKmer(v)
            })
            .collect();
        childs
    }
    fn parents(&self) -> Vec<VecKmer> {
        let (_, prefix) = self.0.split_last().unwrap();
        let parents = BASES
            .iter()
            .map(|first_base| {
                let mut v = prefix.to_vec();
                v.insert(0, *first_base);
                VecKmer(v)
            })
            .collect();
        parents
    }
    fn siblings(&self) -> Vec<VecKmer> {
        let (_, suffix) = self.0.split_first().unwrap();
        BASES
            .iter()
            .map(|first_base| {
                let mut v = suffix.to_vec();
                v.insert(0, *first_base);
                VecKmer(v)
            })
            .collect()
    }
    /// return k+1mer {ACGT}<Kmer> and <Kmer>{ACGT} vector
    fn neighbors(&self) -> Vec<VecKmer> {
        let neighbors: Vec<VecKmer> = VALID_BASES
            .iter()
            .map(|&first_base| {
                let mut v = Vec::new();
                v.push(first_base);
                v.extend_from_slice(&self.0);
                VecKmer(v)
            })
            .chain(VALID_BASES.iter().map(|&last_base| {
                let mut v = Vec::new();
                v.extend_from_slice(&self.0);
                v.push(last_base);
                VecKmer(v)
            }))
            .collect();
        neighbors
    }
    fn extend_first(&self, first_base: u8) -> VecKmer {
        let mut v = Vec::new();
        v.push(first_base);
        v.extend_from_slice(&self.0);
        VecKmer(v)
    }
    fn extend_last(&self, last_base: u8) -> VecKmer {
        let mut v = Vec::new();
        v.extend_from_slice(&self.0);
        v.push(last_base);
        VecKmer(v)
    }
    fn join(&self, other: &VecKmer) -> VecKmer {
        if self.adjacent(other) {
            // self --> other
            self.extend_last(other.last())
        } else if other.adjacent(self) {
            // other --> self
            other.extend_last(self.last())
        } else {
            panic!("cannot join")
        }
    }
    /// check if NNNNNX
    fn is_head(&self) -> bool {
        let k = self.0.len();
        self.0[..k - 1].iter().all(|&x| x == NULL_BASE)
    }
    /// check if XNNNNN
    fn is_tail(&self) -> bool {
        self.0[1..].iter().all(|&x| x == NULL_BASE)
    }
    /// check if not XXXXXN
    fn is_emitable(&self) -> bool {
        *self.0.last().unwrap() != NULL_BASE
    }
    /// check if NXXXXX
    fn is_starting(&self) -> bool {
        self.0[0] == NULL_BASE && self.0[1..].iter().all(|&x| x != NULL_BASE)
    }
    ///
    /// Constructor from slices of u8
    /// (that is slices of DNA bases vector)
    ///
    fn from_bases(bases: &[u8]) -> VecKmer {
        // TODO add base assertion?
        let v = bases.to_vec();
        VecKmer(v)
    }
    ///
    /// Convert back to the raw vector of bases u8
    ///
    fn to_bases(&self) -> Vec<u8> {
        self.0.to_vec()
    }
}

//
// Display
//
impl std::fmt::Display for VecKmer {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // iter returns reference
        for &b in self.0.iter() {
            write!(f, "{}", b as char)?;
        }
        Ok(())
    }
}

//
// Tests
// TODO add test for variable k-mer length
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn veckmer() {
        let a = VecKmer::from_bases(b"ATCGATTAG");
        a.is_emitable();
        println!("{} {}", a, a.len());
    }
}
