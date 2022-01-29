//! TinyKmer definitions
use super::common::{KmerBase, KmerLike};

///
/// Kmer for small k <= 32
/// It can store without heap-allocations
///
#[derive(PartialEq, PartialOrd, Eq, Hash, Clone)]
struct TinyKmer<const K: usize> {
    /// For example,
    ///   s = [s[0]=A, s[1]=C, s[2]=G]
    /// then
    ///   bases = s[0]*(4^2) + s[1]*(4^1) + s[2]*(4^0)
    bases: u64,
    kinds: u64,
}

fn encode_base(base: u8) -> (u64, u64) {
    match base {
        b'A' | b'a' => (0, 0),
        b'C' | b'c' => (1, 0),
        b'G' | b'g' => (2, 0),
        b'T' | b't' => (3, 0),
        _ => (0, 1),
    }
}

fn decode_base(code: u64, kind: u64) -> u8 {
    match kind {
        0 => match code {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => panic!(),
        },
        1 => b'N',
        _ => panic!(),
    }
}

impl TinyKmer<0> {
    fn zero() -> TinyKmer<0> {
        TinyKmer::new(0)
    }
}

impl<const K: usize> TinyKmer<K> {
    fn new(bases: u64) -> TinyKmer<K> {
        TinyKmer { bases, kinds: 0 }
    }
    fn from(bases: &[u8]) -> TinyKmer<K> {
        assert_eq!(bases.len(), K);
        assert!(K <= 32);
        let code = bases
            .iter()
            .map(|&base| encode_base(base).0)
            .fold(0, |acc, code| acc * 4 + code);
        TinyKmer::new(code)
    }
    fn to_vec(&self) -> Vec<u8> {
        let bases = self.bases;
        (0..K)
            .map(|i| {
                let code = (bases >> (2 * i)) % 4;
                decode_base(code, 0)
            })
            .rev()
            .collect()
    }
}

impl<const K: usize> TinyKmer<K>
where
    [(); K + 1]: ,
{
    ///
    /// add base in the first
    /// (BCD, A) -> ABCD
    ///
    fn prepend(&self, base: u8) -> TinyKmer<{ K + 1 }> {
        let (code, kind) = encode_base(base);
        TinyKmer {
            bases: self.bases + (code << (2 * K)),
            kinds: self.kinds + (kind << (2 * K)),
        }
    }
    ///
    /// add base in the last
    /// (ABC, D) -> ABCD
    ///
    fn append(&self, base: u8) -> TinyKmer<{ K + 1 }> {
        let (code, kind) = encode_base(base);
        TinyKmer {
            bases: (self.bases << 2) + code,
            kinds: (self.kinds << 2) + kind,
        }
    }
}

impl<const K: usize> TinyKmer<K>
where
    [(); K - 1]: ,
{
    /// ABCD -> (ABC, D)
    fn pop_last(&self) -> (TinyKmer<{ K - 1 }>, u8) {
        let code = self.bases & 0b11u64;
        let kind = self.kinds & 0b11u64;
        let kmer = TinyKmer {
            bases: self.bases >> 2,
            kinds: self.kinds >> 2,
        };
        let base = decode_base(code, kind);
        (kmer, base)
    }
    /// ABCD -> (BCD, A)
    fn pop_first(&self) -> (TinyKmer<{ K - 1 }>, u8) {
        let code = self.bases & (0b11u64 << (2 * K - 2));
        let kind = self.kinds & (0b11u64 << (2 * K - 2));
        let kmer = TinyKmer {
            bases: self.bases - code,
            kinds: self.kinds - kind,
        };
        let code = code >> (2 * K - 2);
        let kind = kind >> (2 * K - 2);
        let base = decode_base(code, kind);
        (kmer, base)
    }
}

impl<const K: usize> KmerLike for TinyKmer<K>
where
    [(); K - 1]: ,
    [(); K + 1]: ,
{
    type Kp1mer = TinyKmer<{ K + 1 }>;
    type Km1mer = TinyKmer<{ K - 1 }>;
    fn len(&self) -> usize {
        K
    }
    fn first(&self) -> u8 {
        decode_base((self.bases >> (2 * (K - 1))) % 4, 0)
    }
    fn last(&self) -> u8 {
        decode_base(self.bases % 4, 0)
    }
    fn prefix(&self) -> TinyKmer<{ K - 1 }> {
        TinyKmer::new(self.bases >> 2)
    }
    fn suffix(&self) -> TinyKmer<{ K - 1 }> {
        TinyKmer::new(self.bases & !(3 << (2 * K)))
    }
    fn adjacent(&self, other: &TinyKmer<K>) -> bool
    where
        [(); K - 1]: ,
    {
        // for the detail of this bound, see
        // https://github.com/rust-lang/rust/issues/76560
        // this uses nightly feature generic_const_exprs
        self.suffix() == other.prefix()
    }
    fn childs(&self) -> Vec<TinyKmer<K>> {
        unimplemented!();
    }
    fn parents(&self) -> Vec<TinyKmer<K>> {
        unimplemented!();
    }
    fn neighbors(&self) -> Vec<TinyKmer<K>> {
        unimplemented!();
    }
    fn preds(&self) -> Vec<TinyKmer<{ K + 1 }>> {
        unimplemented!();
    }
    fn succs(&self) -> Vec<TinyKmer<{ K + 1 }>> {
        unimplemented!();
    }
    fn join(&self, other: &Self) -> TinyKmer<{ K + 1 }> {
        unimplemented!();
    }
    fn is_head(&self) -> bool {
        unimplemented!();
    }
    fn is_tail(&self) -> bool {
        unimplemented!();
    }
    fn is_emitable(&self) -> bool {
        unimplemented!();
    }
    fn is_starting(&self) -> bool {
        unimplemented!();
    }
    fn extend_first(&self, first_base: u8) -> Self::Kp1mer {
        unimplemented!();
    }
    fn extend_last(&self, last_base: u8) -> Self::Kp1mer {
        unimplemented!();
    }
}

impl<const K: usize> std::fmt::Debug for TinyKmer<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for &b in self.to_vec().iter() {
            write!(f, "{}", b as char)?;
        }
        write!(f, "\nbases={:0>1$b}", self.bases, 64)?;
        write!(f, "\nkinds={:0>1$b}", self.kinds, 64)?;
        Ok(())
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
    fn tinykmer_construct() {
        let m1 = TinyKmer::zero(); // ''
        let m2 = m1.prepend(b'C'); // 'C'
        let m3 = m2.prepend(b'T'); // 'TC'
        let m4 = m3.append(b'G'); // 'TCG'
        println!("{:?}", m1);
        println!("{:?}", m2);
        println!("{:?}", m3);
        println!("{:?}", m4);
        let (m5, b) = m4.pop_first();
        println!("{:?} {}", m5, b as char);
        let (m6, b) = m4.pop_last();
        println!("{:?} {}", m6, b as char);
    }

    #[test]
    fn tinykmer_normal() {
        let va = b"AAAA".to_vec();
        let vb = b"ATCG".to_vec();
        let vc = b"GTACGTA".to_vec();
        let vd = b"TACGTAC".to_vec();
        let a: TinyKmer<4> = TinyKmer::from(&va);
        let b: TinyKmer<4> = TinyKmer::from(&vb);
        let c: TinyKmer<7> = TinyKmer::from(&vc);
        let d: TinyKmer<7> = TinyKmer::from(&vd);

        // convert back
        assert_eq!(a.to_vec(), va);
        assert_eq!(b.to_vec(), vb);
        assert_eq!(c.to_vec(), vc);
        assert_eq!(d.to_vec(), vd);

        // parts
        assert_eq!(b.first(), b'A');
        assert_eq!(b.last(), b'G');
        assert_eq!(b.prefix(), TinyKmer::<3>::from(b"ATC"));
        assert_eq!(b.suffix(), TinyKmer::<3>::from(b"TCG"));
        assert_eq!(b.suffix().to_vec(), b"TCG".to_vec());

        println!("{} {} {} {}", c, d, c.suffix(), d.prefix());
        println!("{:?} {:?}", c.suffix(), d.prefix());
        println!("{} {} {} {}", c, d, c.adjacent(&d), d.adjacent(&c));
        // assert_eq!(c.adjacent(&d), true);
    }

    #[test]
    fn tinykmer_with_n() {
        let ve = b"GTACNN".to_vec();
    }
}
