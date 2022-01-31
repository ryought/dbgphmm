//! VecKmer definitions
use super::common::{KmerBase, KmerLike, NullableKmer};

///
/// Kmer for any k, by using Vec<u8> as a internal struct
///
#[derive(Debug, PartialEq, PartialOrd, Eq, Hash, Clone)]
pub struct VecKmer(pub Vec<u8>);

impl VecKmer {
    ///
    /// Constructor from slices of u8
    /// (that is slices of DNA bases vector)
    ///
    pub fn from(bases: &[u8]) -> VecKmer {
        // TODO add base assertion?
        let v = bases.to_vec();
        VecKmer(v)
    }
    ///
    /// Convert back to the raw vector of bases u8
    ///
    pub fn to_vec(&self) -> Vec<u8> {
        self.0.to_vec()
    }
}

impl NullableKmer for VecKmer {
    /// check if NNNNNN
    fn is_null(&self) -> bool {
        self.0.iter().all(|&x| x == b'N')
    }
}

impl KmerLike for VecKmer {
    type Kp1mer = VecKmer;
    type Km1mer = VecKmer;
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
        let childs = [b'A', b'C', b'G', b'T', b'N']
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
        let parents = [b'A', b'C', b'G', b'T', b'N']
            .iter()
            .map(|first_base| {
                let mut v = prefix.to_vec();
                v.insert(0, *first_base);
                VecKmer(v)
            })
            .collect();
        parents
    }
    /// return k+1mer {ACGT}<Kmer> and <Kmer>{ACGT} vector
    fn neighbors(&self) -> Vec<VecKmer> {
        let bases = [b'A', b'C', b'G', b'T'];
        let neighbors: Vec<VecKmer> = bases
            .iter()
            .map(|&first_base| {
                let mut v = Vec::new();
                v.push(first_base);
                v.extend_from_slice(&self.0);
                VecKmer(v)
            })
            .chain(bases.iter().map(|&last_base| {
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
        self.0[..k - 1].iter().all(|&x| x == b'N')
    }
    /// check if XNNNNN
    fn is_tail(&self) -> bool {
        self.0[1..].iter().all(|&x| x == b'N')
    }
    /// check if not XXXXXN
    fn is_emitable(&self) -> bool {
        *self.0.last().unwrap() != b'N'
    }
    /// check if NXXXXX
    fn is_starting(&self) -> bool {
        self.0[0] == b'N' && self.0[1..].iter().all(|&x| x != b'N')
    }
}

impl KmerBase for VecKmer {
    fn k(&self) -> usize {
        self.0.len()
    }
}

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
        let a = VecKmer::from(b"ATCGATTAG");
        a.is_emitable();
        println!("{} {}", a, a.len());
    }
}
