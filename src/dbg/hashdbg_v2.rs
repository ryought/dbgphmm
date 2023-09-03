//!
//! HashDbg
//!
use crate::common::{CopyNum, Reads, Seq, StyledSequence};
use crate::kmer::kmer::{
    linear_fragment_sequence_to_kmers, sequence_to_kmers, styled_sequence_to_kmers, Kmer, KmerLike,
};
use fnv::{FnvHashMap as HashMap, FnvHashSet as HashSet};
use std::iter::Iterator;

///
/// De Bruijn graph structure that is implemented with a HashMap storing
/// `KmerLike -> CopyNum` mapping.
///
#[derive(Debug)]
pub struct HashDbg<K: KmerLike> {
    k: usize,
    kmers: HashMap<K, CopyNum>,
}

///
/// Basic operations
///
impl<K: KmerLike> HashDbg<K> {
    /// Create new HashDbg with no k-mers
    pub fn new(k: usize) -> HashDbg<K> {
        HashDbg {
            k,
            kmers: HashMap::default(),
        }
    }
    /// size of k-mer in HashDbg
    pub fn k(&self) -> usize {
        self.k
    }
    ///
    pub fn set(&mut self, kmer: K, copy_num: CopyNum) {
        assert_eq!(kmer.k(), self.k());
        if copy_num > 0 {
            self.kmers.insert(kmer, copy_num);
        } else {
            self.kmers.remove(&kmer);
        }
    }
    ///
    pub fn get(&self, kmer: &K) -> CopyNum {
        assert_eq!(kmer.k(), self.k());
        match self.kmers.get(&kmer) {
            Some(&copy_num) => copy_num,
            _ => 0,
        }
    }
    /// Add k-mer count
    pub fn add(&mut self, kmer: K, copy_num: CopyNum) {
        let copy_num_old = self.get(&kmer);
        self.set(kmer, copy_num + copy_num_old);
    }
    /// Check k-mer (edge) exists in DBG
    pub fn has(&self, kmer: &K) -> bool {
        assert_eq!(kmer.k(), self.k());
        self.get(kmer) > 0
    }
    /// Get a list of edges (k-mers)
    pub fn edges(&self) -> Vec<K> {
        self.kmers.keys().cloned().collect()
    }
    /// Get a list of nodes (k-1-mers)
    pub fn nodes(&self) -> Vec<K> {
        let mut nodes = HashSet::default();
        for kmer in self.kmers.keys() {
            nodes.insert(kmer.prefix());
            nodes.insert(kmer.suffix());
        }
        nodes.into_iter().collect()
    }
    pub fn kmers(&self) -> impl Iterator<Item = &K> {
        self.kmers.keys()
    }
    /// Source node (k-1-mer) from edge (k-mer)
    pub fn source(&self, kmer: &K) -> K {
        assert!(self.has(kmer));
        kmer.prefix()
    }
    /// Target node (k-1-mer) from edge (k-mer)
    pub fn target(&self, kmer: &K) -> K {
        assert!(self.has(kmer));
        kmer.suffix()
    }
    /// List of incoming edges of node (k-1-mer)
    /// km1mer XX -> kmers [YXX] whose suffix is the km1mer
    pub fn edges_in(&self, km1mer: &K) -> Vec<K> {
        assert_eq!(km1mer.k(), self.k() - 1);
        km1mer
            .preds()
            .into_iter()
            .filter(|kmer| self.has(kmer))
            .collect()
    }
    /// List of outgoing edges of node (k-1-mer)
    /// km1mer XX -> kmers [XXY] whose prefix is the km1mer
    pub fn edges_out(&self, km1mer: &K) -> Vec<K> {
        assert_eq!(km1mer.k(), self.k() - 1);
        km1mer
            .succs()
            .into_iter()
            .filter(|kmer| self.has(kmer))
            .collect()
    }
    /// Deprecated
    pub fn childs(&self, kmer: &K) -> Vec<K> {
        assert_eq!(kmer.k(), self.k());
        kmer.childs()
            .into_iter()
            .filter(|child| self.has(child))
            .collect()
    }
    /// Deprecated
    pub fn parents(&self, kmer: &K) -> Vec<K> {
        assert_eq!(kmer.k(), self.k());
        kmer.parents()
            .into_iter()
            .filter(|parent| self.has(parent))
            .collect()
    }
    /// Deprecated
    pub fn siblings(&self, kmer: &K) -> Vec<K> {
        assert_eq!(kmer.k(), self.k());
        kmer.siblings()
            .into_iter()
            .filter(|sibling| self.has(sibling))
            .collect()
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
    /// add all kmers in linear fragment seq (without NNN kmers)
    ///
    pub fn add_fragment_seq(&mut self, seq: &[u8]) {
        for kmer in linear_fragment_sequence_to_kmers(seq, self.k()) {
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
            assert!(!d.has(kmer));
            d.set(kmer.clone(), *copy_num);
        }
        d
    }
    pub fn from_seq<S: Seq>(k: usize, seq: &S) -> Self {
        let mut d = HashDbg::new(k);
        d.add_seq(seq.as_ref());
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
    /// Construct HashDbg from seqs regarding as Fragments
    /// Ignoring heads and tails (like nnnA and Gnnn).
    ///
    pub fn from_fragment_seqs<T>(k: usize, seqs: T) -> Self
    where
        T: IntoIterator,
        T::Item: Seq,
    {
        let mut d = HashDbg::new(k);
        for seq in seqs {
            let seq = seq.as_ref();
            if seq.len() >= k {
                d.add_fragment_seq(seq);
            }
        }
        d
    }
    pub fn from_styled_seqs<T>(k: usize, seqs: T) -> Self
    where
        T: IntoIterator,
        T::Item: AsRef<StyledSequence>,
    {
        let mut d = HashDbg::new(k);
        for seq in seqs {
            d.add_styled_sequence(seq.as_ref());
        }
        d
    }
}

impl<K: KmerLike> HashDbg<K> {
    ///
    /// to kmer count profile
    ///
    pub fn to_kmer_profile(&self) -> HashMap<K, CopyNum> {
        self.kmers.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kmer::veckmer::{kmer, VecKmer};

    #[test]
    fn hashdbg_v2_new() {
        let mut hd: HashDbg<VecKmer> = HashDbg::new(4);

        // add get
        assert_eq!(hd.get(&kmer(b"ATCG")), 0);
        hd.add(kmer(b"ATCG"), 1);
        hd.add(kmer(b"TTTT"), 2);
        assert_eq!(hd.get(&kmer(b"TTTT")), 2);
        hd.add(kmer(b"TTTT"), 1);
        assert_eq!(hd.get(&kmer(b"TTTT")), 3);
        assert_eq!(hd.get(&kmer(b"ATCG")), 1);
        println!("{:?}", hd);

        // is_exists
        assert_eq!(hd.has(&kmer(b"ATCG")), true);
        assert_eq!(hd.has(&kmer(b"ATCT")), false);

        // kmers
        let kmers: Vec<VecKmer> = hd.kmers().cloned().collect();
        assert_eq!(kmers, vec![kmer(b"TTTT"), kmer(b"ATCG")]);
        assert_eq!(hd.edges(), vec![kmer(b"TTTT"), kmer(b"ATCG")]);
        assert_eq!(hd.nodes(), vec![kmer(b"ATC"), kmer(b"TCG"), kmer(b"TTT")]);
        assert_eq!(hd.edges_out(&kmer(b"ATC")), vec![kmer(b"ATCG")]);
        assert_eq!(hd.edges_in(&kmer(b"ATC")), vec![]);
        assert_eq!(hd.edges_out(&kmer(b"TTT")), vec![kmer(b"TTTT")]);
        assert_eq!(hd.edges_in(&kmer(b"TTT")), vec![kmer(b"TTTT")]);

        // set and delete
        hd.set(Kmer::from_bases(b"TTTT"), 0);
        assert_eq!(hd.has(&Kmer::from_bases(b"TTTT")), false);

        assert_eq!(hd.edges(), vec![kmer(b"ATCG")]);
        assert_eq!(hd.nodes(), vec![kmer(b"ATC"), kmer(b"TCG")]);
    }

    #[test]
    fn hashdbg_v2_seq() {
        let mut hd: HashDbg<VecKmer> = HashDbg::from_seq(4, b"ATCGATTCGAT");
        println!("{}", hd);

        // childs parents siblings
        assert_eq!(hd.childs(&kmer(b"ATCG")), vec![kmer(b"TCGA")]);
        assert_eq!(hd.childs(&kmer(b"ATCC")), vec![]);
        assert_eq!(hd.parents(&kmer(b"CGAT")), vec![kmer(b"TCGA")]);
        assert_eq!(
            hd.childs(&kmer(b"CGAT")),
            vec![kmer(b"GATT"), kmer(b"GATn")]
        );
        assert_eq!(hd.childs(&kmer(b"Tnnn")), vec![kmer(b"nnnA")]);
        assert_eq!(
            hd.siblings(&kmer(b"ATCG")),
            vec![kmer(b"ATCG"), kmer(b"TTCG")]
        );
        assert_eq!(hd.siblings(&kmer(b"nnnA")), vec![kmer(b"nnnA")]);

        assert_eq!(hd.edges().len(), 12);
        assert_eq!(hd.nodes().len(), 11);

        for e in hd.edges_out(&kmer(b"nnn")) {
            println!("out {}", e);
        }
        assert_eq!(hd.edges_out(&kmer(b"nnn")), vec![kmer(b"nnnA")]);
        assert_eq!(hd.edges_in(&kmer(b"nnn")), vec![kmer(b"Tnnn")]);
        assert_eq!(hd.edges_out(&kmer(b"TCG")), vec![kmer(b"TCGA")]);
        assert_eq!(
            hd.edges_out(&kmer(b"GAT")),
            vec![kmer(b"GATT"), kmer(b"GATn")]
        );
        assert_eq!(hd.get(&kmer(b"TCGA")), 2);
        assert_eq!(hd.source(&kmer(b"TCGA")), kmer(b"TCG"));
        assert_eq!(hd.target(&kmer(b"TCGA")), kmer(b"CGA"));
    }
    #[test]
    fn hashdbg_v2_profile() {
        let hd: HashDbg<VecKmer> = HashDbg::from_seq(4, b"ATCGATTCGAT");
        let p = hd.to_kmer_profile();
        for (k, v) in p.iter() {
            println!("{} {}", k, v);
        }

        // answer
        let mut p2 = HashMap::default();
        p2.insert(kmer(b"ATTC"), 1);
        p2.insert(kmer(b"TTCG"), 1);
        p2.insert(kmer(b"nnAT"), 1);
        p2.insert(kmer(b"CGAT"), 2);
        p2.insert(kmer(b"TCGA"), 2);
        p2.insert(kmer(b"nATC"), 1);
        p2.insert(kmer(b"GATT"), 1);
        p2.insert(kmer(b"GATn"), 1);
        p2.insert(kmer(b"nnnA"), 1);
        p2.insert(kmer(b"ATnn"), 1);
        p2.insert(kmer(b"ATCG"), 1);
        p2.insert(kmer(b"Tnnn"), 1);

        assert_eq!(p, p2);
    }
}
