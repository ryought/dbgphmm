//! VecKmer definitions
use super::common::{KmerBase, KmerLike};

///
/// Kmer for any k
///
#[derive(Debug, PartialEq, PartialOrd, Eq, Hash, Clone)]
struct VecKmer(Vec<u8>);

impl VecKmer {
    fn from(bases: &[u8]) -> VecKmer {
        let v = bases.to_vec();
        VecKmer(v)
    }
}

impl KmerBase for VecKmer {
    type Kp1mer = VecKmer;
    type Km1mer = VecKmer;
    fn k(&self) -> usize {
        self.0.len()
    }
    fn prefix(&self) -> VecKmer {
        let (_, prefix) = self.0.split_last().unwrap();
        VecKmer(prefix.to_vec())
    }
    fn suffix(&self) -> VecKmer {
        let (_, suffix) = self.0.split_first().unwrap();
        VecKmer(suffix.to_vec())
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
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn veckmer() {
        let a = VecKmer::from(b"ATCGATTAG");
        println!("{} {}", a, a.k());
    }
}
