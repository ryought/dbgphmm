//!
//! Kmer definitions
//!

// traits
trait Kmer: std::marker::Sized {
    // basic
    fn k(&self) -> usize;
    /*
    // construction
    fn from(bases: &[u8]) -> Self;
    fn to_vec(&self) -> Vec<u8>;
    // fn to_veckmer(&self) -> VecKmer;
    // fn to_tinykmer(&self) -> TinyKmer;
    // adjacency
    fn is_adjacent<T: Kmer>(&self, other: &T) -> bool;
    // parts
    fn first(&self) -> u8;
    fn last(&self) -> u8;
    fn prefix(&self) -> Self;
    fn suffix(&self) -> Self;
    fn childs(&self) -> Vec<Self>;
    fn parents(&self) -> Vec<Self>;
    */
}

///
/// Kmer for k > 32
///
#[derive(Debug, PartialEq, PartialOrd, Eq, Hash, Clone)]
struct VecKmer(Vec<u8>);

impl VecKmer {
    fn from(bases: &[u8]) -> VecKmer {
        let v = bases.to_vec();
        VecKmer(v)
    }
}

impl Kmer for VecKmer {
    fn k(&self) -> usize {
        self.0.len()
    }
}

///
/// Kmer for k <= 32
/// It can store without heap-allocations
///
#[derive(Debug, PartialEq, PartialOrd, Eq, Hash, Clone)]
struct TinyKmer<const K: usize> {
    bases: u64,
}

impl<const K: usize> TinyKmer<K> {
    fn new() -> TinyKmer<K> {
        TinyKmer { bases: 1 }
    }
    fn from(bases: &[u8]) -> TinyKmer<K> {
        assert_eq!(bases.len(), K);
        TinyKmer { bases: 10 }
    }
    fn prefix(&self) -> TinyKmer<{ K - 1 }> {
        TinyKmer::<{ K - 1 }> { bases: 0 }
    }
}

impl<const K: usize> Kmer for TinyKmer<K> {
    fn k(&self) -> usize {
        K
    }
}

//
// Tests
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tinykmer() {
        let kmer1: TinyKmer<10> = TinyKmer::new();
        println!("{:?}", kmer1);
        println!("{}", kmer1.k());

        let kmer2 = kmer1.prefix();
        println!("{:?}", kmer2);
        println!("{}", kmer2.k());
    }
}
