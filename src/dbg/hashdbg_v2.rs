//!
//! HashDbg
//!
use crate::common::{CopyNum, Reads, Seq, StyledSequence};
use crate::kmer::kmer::{sequence_to_kmers, styled_sequence_to_kmers, Kmer, KmerLike};
use fnv::FnvHashMap as HashMap;
use std::iter::Iterator;

///
/// De Bruijn graph structure that is implemented with a HashMap storing
/// `KmerLike -> CopyNum` mapping.
///
#[derive(Debug)]
pub struct HashDbg<K: KmerLike> {
    k: usize,
    store: HashMap<K, CopyNum>,
}

///
/// Basic operations
///
impl<K: KmerLike> HashDbg<K> {
    pub fn new(k: usize) -> HashDbg<K> {
        HashDbg {
            k,
            store: HashMap::default(),
        }
    }
    pub fn k(&self) -> usize {
        self.k
    }
    pub fn set(&mut self, kmer: K, copy_num: CopyNum) {
        assert_eq!(kmer.k(), self.k());
        if copy_num > 0 {
            self.store.insert(kmer, copy_num);
        } else {
            self.store.remove(&kmer);
        }
    }
    pub fn get(&self, kmer: &K) -> CopyNum {
        assert_eq!(kmer.k(), self.k());
        match self.store.get(&kmer) {
            Some(&copy_num) => copy_num,
            _ => 0,
        }
    }
    pub fn add(&mut self, kmer: K, copy_num: CopyNum) {
        let copy_num_old = self.get(&kmer);
        self.set(kmer, copy_num + copy_num_old);
    }
    pub fn is_exists(&self, kmer: &K) -> bool {
        assert_eq!(kmer.k(), self.k());
        self.get(kmer) > 0
    }
    pub fn kmers(&self) -> impl Iterator<Item = &K> {
        self.store.keys()
    }
    pub fn childs(&self, kmer: &K) -> Vec<K> {
        assert_eq!(kmer.k(), self.k());
        kmer.childs()
            .into_iter()
            .filter(|child| self.is_exists(child))
            .collect()
    }
    pub fn parents(&self, kmer: &K) -> Vec<K> {
        assert_eq!(kmer.k(), self.k());
        kmer.parents()
            .into_iter()
            .filter(|parent| self.is_exists(parent))
            .collect()
    }
    pub fn siblings(&self, kmer: &K) -> Vec<K> {
        assert_eq!(kmer.k(), self.k());
        kmer.siblings()
            .into_iter()
            .filter(|sibling| self.is_exists(sibling))
            .collect()
    }
    pub fn is_consistent(&self) -> bool {
        self.kmers().all(|kmer| {
            // (sum of child kmers) =? (sum of sibling kmers)
            let n1: CopyNum = self.childs(&kmer).iter().map(|child| self.get(child)).sum();
            let n2: CopyNum = self
                .siblings(&kmer)
                .iter()
                .map(|sibling| self.get(sibling))
                .sum();
            n1 == n2
        })
    }
    ///
    /// add all kmers in linear seq (with leading/trailing NNN kmers)
    ///
    pub fn add_seq(&mut self, seq: &[u8]) {
        for kmer in sequence_to_kmers(seq, self.k()) {
            self.add(kmer, 1);
        }
    }
    ///
    /// add all kmers in styled sequence
    ///
    pub fn add_styled_sequence(&mut self, s: &StyledSequence) {
        for kmer in styled_sequence_to_kmers(s, self.k()) {
            self.add(kmer, 1);
        }
    }
}

//
// Display
//
impl<K: KmerLike + std::fmt::Display> std::fmt::Display for HashDbg<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // iter returns reference
        for kmer in self.kmers() {
            writeln!(f, "{} {}", kmer, self.get(&kmer))?;
        }
        Ok(())
    }
}

///
/// Constructors
///
impl<K: KmerLike> HashDbg<K> {
    pub fn from_profile(k: usize, profile: &[(K, CopyNum)]) -> Self {
        let mut d = HashDbg::new(k);
        for (kmer, copy_num) in profile.iter() {
            assert!(!d.is_exists(kmer));
            d.set(kmer.clone(), *copy_num);
        }
        d
    }
    pub fn from_seq(k: usize, seq: &[u8]) -> Self {
        let mut d = HashDbg::new(k);
        d.add_seq(seq);
        d
    }
    pub fn from_seqs<T>(k: usize, seqs: T) -> Self
    where
        T: IntoIterator,
        T::Item: Seq,
    {
        let mut d = HashDbg::new(k);
        for seq in seqs {
            let seq = seq.as_ref();
            // ignore read if it is shorter than k
            if seq.len() >= k {
                d.add_seq(seq);
            }
        }
        d
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::veckmer::VecKmer;

    #[test]
    fn hashdbg_v2_new() {
        let mut hd: HashDbg<VecKmer> = HashDbg::new(4);

        // add get
        assert_eq!(hd.get(&Kmer::from_bases(b"ATCG")), 0);
        hd.add(Kmer::from_bases(b"ATCG"), 1);
        hd.add(Kmer::from_bases(b"TTTT"), 2);
        assert_eq!(hd.get(&Kmer::from_bases(b"TTTT")), 2);
        hd.add(Kmer::from_bases(b"TTTT"), 1);
        assert_eq!(hd.get(&Kmer::from_bases(b"TTTT")), 3);
        assert_eq!(hd.get(&Kmer::from_bases(b"ATCG")), 1);
        println!("{:?}", hd);

        // is_exists
        assert_eq!(hd.is_exists(&Kmer::from_bases(b"ATCG")), true);
        assert_eq!(hd.is_exists(&Kmer::from_bases(b"ATCT")), false);

        // kmers
        let kmers: Vec<VecKmer> = hd.kmers().cloned().collect();
        assert_eq!(
            kmers,
            vec![Kmer::from_bases(b"TTTT"), Kmer::from_bases(b"ATCG"),]
        );

        // set and delete
        hd.set(Kmer::from_bases(b"TTTT"), 0);
        assert_eq!(hd.is_exists(&Kmer::from_bases(b"TTTT")), false);

        // not consistent
        assert!(!hd.is_consistent());
    }

    #[test]
    fn hashdbg_v2_seq() {
        let mut hd: HashDbg<VecKmer> = HashDbg::from_seq(4, b"ATCGATTCGAT");
        println!("{}", hd);

        // childs parents siblings
        assert_eq!(
            hd.childs(&Kmer::from_bases(b"ATCG")),
            vec![Kmer::from_bases(b"TCGA")]
        );
        assert_eq!(hd.childs(&Kmer::from_bases(b"ATCC")), vec![]);
        assert_eq!(
            hd.parents(&Kmer::from_bases(b"CGAT")),
            vec![Kmer::from_bases(b"TCGA")]
        );
        assert_eq!(
            hd.childs(&Kmer::from_bases(b"CGAT")),
            vec![Kmer::from_bases(b"GATT"), Kmer::from_bases(b"GATn")]
        );
        assert_eq!(
            hd.childs(&Kmer::from_bases(b"Tnnn")),
            vec![Kmer::from_bases(b"nnnA")]
        );
        assert_eq!(
            hd.siblings(&Kmer::from_bases(b"ATCG")),
            vec![Kmer::from_bases(b"ATCG"), Kmer::from_bases(b"TTCG")]
        );
        assert_eq!(
            hd.siblings(&Kmer::from_bases(b"nnnA")),
            vec![Kmer::from_bases(b"nnnA")]
        );

        // consistency
        assert!(hd.is_consistent());
    }
}
