//!
//! Kmer definitions
//!
use super::kmer::Kmer;

pub trait KmerLike: std::marker::Sized {
    fn len(&self) -> usize;
    fn k(&self) -> usize {
        self.len()
    }
    fn adjacent(&self, other: &Kmer) -> bool;
    fn first(&self) -> u8;
    fn last(&self) -> u8;
    fn prefix(&self) -> Kmer;
    fn suffix(&self) -> Kmer;
    fn childs(&self) -> Vec<Kmer>;
    fn parents(&self) -> Vec<Kmer>;
    fn neighbors(&self) -> Vec<Kmer>;
    fn preds(&self) -> Vec<Kmer>;
    fn succs(&self) -> Vec<Kmer>;
    fn extend_first(&self, first_base: u8) -> Kmer;
    fn extend_last(&self, last_base: u8) -> Kmer;
    fn join(&self, other: &Kmer) -> Kmer;
    fn is_head(&self) -> bool;
    fn is_tail(&self) -> bool;
    fn is_emitable(&self) -> bool;
    fn is_starting(&self) -> bool;
    // construction
    // fn from(bases: &[u8]) -> Self;
    // fn to_vec(&self) -> Vec<u8>;
}

trait KmerV2 {
    fn k(&self) -> usize;
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

impl KmerV2 for VecKmer {
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
    fn head(&self) -> u8 {
        decode_base((self.bases >> (2 * (K - 1))) % 4)
    }
    fn tail(&self) -> u8 {
        decode_base(self.bases % 4)
    }
    fn prefix(&self) -> TinyKmer<{ K - 1 }> {
        TinyKmer {
            bases: self.bases >> 2,
        }
    }
    fn suffix(&self) -> TinyKmer<{ K - 1 }> {
        TinyKmer {
            bases: self.bases & !(3 << (2 * K)),
        }
    }
    fn is_adjacent(&self, other: &TinyKmer<K>) -> bool
    where
        [(); K - 1]: ,
    {
        // for the detail of this bound, see
        // https://github.com/rust-lang/rust/issues/76560
        // this uses nightly feature generic_const_exprs
        self.suffix() == other.prefix()
    }
}

impl<const K: usize> KmerV2 for TinyKmer<K> {
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
        println!("{:?}", vb);
        assert_eq!(b.head(), b'A');
        assert_eq!(b.tail(), b'G');
        println!("{}", b.prefix());
        println!("{}", b.suffix());
    }
}
