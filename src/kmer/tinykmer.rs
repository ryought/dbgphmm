//! TinyKmer definitions
use super::common::{KmerBase, KmerLike};

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
    fn _prefix(&self) -> TinyKmer<{ K - 1 }> {
        TinyKmer {
            bases: self.bases >> 2,
        }
    }
    fn _suffix(&self) -> TinyKmer<{ K - 1 }> {
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
        self._suffix() == other._prefix()
    }
}

impl<const K: usize> KmerBase for TinyKmer<K>
where
    [(); K - 1]: ,
    [(); K + 1]: ,
{
    type Kp1mer = TinyKmer<{ K + 1 }>;
    type Km1mer = TinyKmer<{ K - 1 }>;
    fn k(&self) -> usize {
        K
    }
    fn prefix(&self) -> TinyKmer<{ K - 1 }> {
        self._prefix()
    }
    fn suffix(&self) -> TinyKmer<{ K - 1 }> {
        self._suffix()
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
