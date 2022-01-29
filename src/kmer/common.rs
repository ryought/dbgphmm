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

impl Kmer for VecKmer {
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

///
/// Kmer for small k <= 32
/// It can store without heap-allocations
///
#[derive(Debug, PartialEq, PartialOrd, Eq, Hash, Clone)]
struct TinyKmer<const K: usize> {
    bases: u64,
}

fn encode_base(base: u8) -> u64 {
    match base {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => panic!(),
    }
}
fn decode_base(code: u64) -> u8 {
    match code {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => panic!(),
    }
}

impl<const K: usize> TinyKmer<K> {
    fn from(bases: &[u8]) -> TinyKmer<K> {
        assert_eq!(bases.len(), K);
        assert!(K <= 32);
        let code = bases
            .iter()
            .map(|&base| encode_base(base))
            .fold(0, |acc, code| acc * 4 + code);
        TinyKmer { bases: code }
    }
    fn to_vec(&self) -> Vec<u8> {
        let bases = self.bases;
        (0..K)
            .map(|i| {
                let code = (bases >> (2 * i)) % 4;
                decode_base(code)
            })
            .rev()
            .collect()
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

impl<const K: usize> std::fmt::Display for TinyKmer<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for &b in self.to_vec().iter() {
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

    #[test]
    fn tinykmer() {
        let va = b"AAAA".to_vec();
        let vb = b"ATCG".to_vec();
        let vc = b"GTACGTA".to_vec();
        let a: TinyKmer<4> = TinyKmer::from(&va);
        let b: TinyKmer<4> = TinyKmer::from(&vb);
        let c: TinyKmer<7> = TinyKmer::from(&vc);
        assert_eq!(a.to_vec(), va);
        assert_eq!(b.to_vec(), vb);
        assert_eq!(c.to_vec(), vc);
    }
}
