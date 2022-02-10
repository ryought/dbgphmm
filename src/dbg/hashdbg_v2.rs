//!
//! HashDbg
//!
use crate::common::CopyNum;
use crate::kmer::kmer::{ending_kmers, starting_kmers, Kmer, KmerLike};
use fnv::FnvHashMap as HashMap;
use std::iter::Iterator;

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
    pub fn is_consistent(&self) -> bool {
        for kmer in self.kmers() {
            // (sum of child kmers) =? (sum of sibling kmers)
        }
        true
    }
    pub fn add_seq(&mut self, seq: &[u8]) {}
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
        // TODO
        d
    }
}
