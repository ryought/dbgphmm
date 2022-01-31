//! TinyKmer definitions
use super::common::{KmerLike, NullableKmer};
use super::quadarray::QuadArray;

///
/// Kmer for small k <= 32
/// It can store without heap-allocations
///
#[derive(PartialEq, PartialOrd, Eq, Hash, Clone, Copy)]
pub struct TinyKmer<const K: usize> {
    pub codes: QuadArray,
    pub kinds: QuadArray,
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
    // constructors
    fn new(bases: u64) -> TinyKmer<K> {
        TinyKmer {
            codes: QuadArray::empty(),
            kinds: QuadArray::empty(),
        }
    }
    fn empty() -> TinyKmer<K> {
        TinyKmer {
            codes: QuadArray::empty(),
            kinds: QuadArray::empty(),
        }
    }
    pub fn from(bases: &[u8]) -> TinyKmer<K> {
        assert_eq!(bases.len(), K);
        assert!(K <= 32);
        let mut kmer = TinyKmer::empty();
        for i in 0..K {
            kmer.set(i, bases[i]);
        }
        kmer
    }
    fn to_vec(&self) -> Vec<u8> {
        let mut bases = vec![0; K];
        for i in 0..K {
            bases[i] = self.get(i);
        }
        bases
    }
    // element-wise functions
    fn get(&self, index: usize) -> u8 {
        let code = self.codes.get(index);
        let kind = self.kinds.get(index);
        decode_base(code, kind)
    }
    fn set(&mut self, index: usize, base: u8) {
        let (code, kind) = encode_base(base);
        self.codes.set(index, code);
        self.kinds.set(index, kind);
    }
    fn shift_back(&mut self) {
        self.codes.shift_back();
        self.kinds.shift_back();
    }
    fn shift_front(&mut self) {
        self.codes.shift_front();
        self.kinds.shift_front();
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
        // copy to new kmer that will be returned
        // with the same content, but different type
        let mut kmer: TinyKmer<{ K + 1 }> = TinyKmer {
            codes: self.codes,
            kinds: self.kinds,
        };
        kmer.shift_back();
        kmer.set(0, base);
        kmer
    }
    ///
    /// add base in the last
    /// (ABC, D) -> ABCD
    ///
    fn append(&self, base: u8) -> TinyKmer<{ K + 1 }> {
        // copy to new kmer that will be returned
        // with the same content, but different type
        let mut kmer: TinyKmer<{ K + 1 }> = TinyKmer {
            codes: self.codes,
            kinds: self.kinds,
        };
        kmer.set(K, base);
        kmer
    }
}

impl<const K: usize> TinyKmer<K>
where
    [(); K - 1]: ,
{
    /// ABCD -> (ABC, D)
    fn pop_last(&self) -> (TinyKmer<{ K - 1 }>, u8) {
        // get the last element
        let code = self.codes.get(K - 1);
        let kind = self.kinds.get(K - 1);
        let base = decode_base(code, kind);
        let mut kmer: TinyKmer<{ K - 1 }> = TinyKmer {
            codes: self.codes,
            kinds: self.kinds,
        };
        // delete the last element
        kmer.codes.set(K - 1, 0);
        kmer.kinds.set(K - 1, 0);
        (kmer, base)
    }
    /// ABCD -> (BCD, A)
    fn pop_first(&self) -> (TinyKmer<{ K - 1 }>, u8) {
        let code = self.codes.get(0);
        let kind = self.kinds.get(0);
        let base = decode_base(code, kind);
        let mut kmer: TinyKmer<{ K - 1 }> = TinyKmer {
            codes: self.codes,
            kinds: self.kinds,
        };
        // move all elements toward front
        kmer.codes.shift_front();
        kmer.kinds.shift_front();
        (kmer, base)
    }
}

impl<const K: usize> NullableKmer for TinyKmer<K> {
    fn is_null(&self) -> bool {
        true
    }
}

impl<const K: usize> KmerLike for TinyKmer<K>
where
    [(); K - 1]: ,
    [(); K + 1]: ,
{
    // for the detail of this bound in where, see
    // https://github.com/rust-lang/rust/issues/76560
    // this uses nightly feature generic_const_exprs
    type Kp1mer = TinyKmer<{ K + 1 }>;
    type Km1mer = TinyKmer<{ K - 1 }>;
    fn len(&self) -> usize {
        K
    }
    fn k(&self) -> usize {
        K
    }
    fn first(&self) -> u8 {
        self.get(0)
    }
    fn last(&self) -> u8 {
        self.get(K - 1)
    }
    fn prefix(&self) -> TinyKmer<{ K - 1 }> {
        let (kmer, _) = self.pop_last();
        kmer
    }
    fn suffix(&self) -> TinyKmer<{ K - 1 }> {
        let (kmer, _) = self.pop_first();
        kmer
    }
    fn adjacent(&self, other: &TinyKmer<K>) -> bool {
        self.suffix() == other.prefix()
    }
    fn childs(&self) -> Vec<TinyKmer<K>> {
        let suffix = self.suffix();
        [b'A', b'C', b'G', b'T', b'N']
            .iter()
            .map(|&last_base| {
                let mut child: TinyKmer<K> = TinyKmer::empty();
                child.codes = suffix.codes;
                child.kinds = suffix.kinds;
                child.set(K - 1, last_base);
                child
            })
            .collect()
    }
    fn parents(&self) -> Vec<TinyKmer<K>> {
        let prefix = self.prefix();
        [b'A', b'C', b'G', b'T', b'N']
            .iter()
            .map(|&first_base| {
                let mut parent: TinyKmer<K> = TinyKmer::empty();
                parent.codes = prefix.codes;
                parent.kinds = prefix.kinds;
                parent.shift_back();
                parent.set(0, first_base);
                parent
            })
            .collect()
    }
    fn extend_first(&self, first_base: u8) -> Self::Kp1mer {
        self.prepend(first_base)
    }
    fn extend_last(&self, last_base: u8) -> Self::Kp1mer {
        self.append(last_base)
    }
}

impl<const K: usize> std::fmt::Debug for TinyKmer<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for &b in self.to_vec().iter() {
            write!(f, "{}", b as char)?;
        }
        write!(f, "\ncodes={:0>1$b}", self.codes.0, 64)?;
        write!(f, "\nkinds={:0>1$b}", self.kinds.0, 64)?;
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
        let (m5, b5) = m4.pop_first();
        let (m6, b6) = m4.pop_last();

        assert_eq!(m1, TinyKmer::from(b""));
        assert_eq!(m2, TinyKmer::from(b"C"));
        assert_eq!(m3, TinyKmer::from(b"TC"));
        assert_eq!(m4, TinyKmer::from(b"TCG"));
        assert_eq!(m5, TinyKmer::from(b"CG"));
        assert_eq!(b5, b'T');
        assert_eq!(m6, TinyKmer::from(b"TC"));
        assert_eq!(b6, b'G');

        let m7 = m4.append(b'N');
        assert_eq!(m7, TinyKmer::from(b"TCGN"));
        let (m8, b8) = m7.pop_last();
        assert_eq!(m8, m4);
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
    }

    #[test]
    fn tinykmer_recursive() {
        let va = b"ATCGT".to_vec();
        let a: TinyKmer<5> = TinyKmer::from(&va);
        let b = a.suffix();
        println!("a={}", a);
        println!("b={}", b);
        println!("{}", b.prepend(b'A'));
    }

    /*
    #[test]
    fn tinykmer_as_kmerlike() {
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
    */
}
