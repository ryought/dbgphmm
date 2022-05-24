//!
//! Comparison methods of two de Bruijn graphs.
//!
use super::dbg::{Dbg, DbgEdge, DbgNode};
use crate::common::{CopyNum, Seq, SeqStyle, Sequence};
use crate::dbg::hashdbg_v2::HashDbg;
use crate::hist::Hist;
use crate::kmer::common::linear_sequence_to_kmers;
use crate::kmer::{KmerLike, NullableKmer};

///
/// Result of comparing two Dbgs.
///
#[derive(Clone, Debug, Default)]
pub struct CompareResult {
    pub n_true: usize,
    pub n_error: usize,
}

///
/// Result of comparing two Dbgs along with a seq.
///
#[derive(Clone, Debug, Default)]
pub struct CompareWithSeqResult<K: KmerLike>(pub Vec<CompareWithSeqResultForBase<K>>);

#[derive(Clone, Debug, Default)]
pub struct CompareWithSeqResultForBase<K: KmerLike> {
    kmer: K,
    copy_num_true: CopyNum,
    copy_num_self: CopyNum,
}

impl<K: KmerLike> CompareWithSeqResultForBase<K> {
    pub fn is_correct(&self) -> bool {
        self.copy_num_true == self.copy_num_self
    }
}

impl<K: KmerLike> std::fmt::Display for CompareWithSeqResult<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for (i, r) in self.0.iter().enumerate() {
            writeln!(
                f,
                "{}\t{}\t{}\t{}",
                i, r.kmer, r.copy_num_true, r.copy_num_self
            )?;
        }
        Ok(())
    }
}

///
/// result struct of `check_kmer_existence_with_seq`
///
#[derive(Clone, Debug)]
pub struct KmerExistenceResult<K: KmerLike> {
    pub n_exists: usize,
    pub n_not_exists: usize,
    pub kmers_not_exists: Vec<K>,
}

impl<K: KmerLike> KmerExistenceResult<K> {
    pub fn new() -> Self {
        KmerExistenceResult {
            n_exists: 0,
            n_not_exists: 0,
            kmers_not_exists: Vec::new(),
        }
    }
}

///
/// Histgram of copy_nums of kmers with the specified copy_num.
///
#[derive(Clone, Debug)]
pub struct KmerHists(Vec<Hist>);

impl KmerHists {
    ///
    /// Create a new kmer hists struct
    ///
    pub fn new(max_copy_num: usize) -> Self {
        let hists = vec![Hist::new(); max_copy_num + 1];
        KmerHists(hists)
    }
    ///
    /// get the maximum copy number
    ///
    pub fn max_copy_num(&self) -> usize {
        self.0.len() - 1
    }
    ///
    ///
    ///
    pub fn get_hist(&self, copy_num: usize) -> &Hist {
        &self.0[copy_num]
    }
    ///
    ///
    ///
    pub fn get_hist_mut(&mut self, copy_num: usize) -> &mut Hist {
        &mut self.0[copy_num]
    }
}

impl std::fmt::Display for KmerHists {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        (0..=self.max_copy_num()).try_for_each(|copy_num| {
            let hist = self.get_hist(copy_num);
            write!(
                f,
                "x{}={}{}",
                copy_num,
                hist,
                if copy_num != self.max_copy_num() {
                    ";"
                } else {
                    ""
                }
            )
        })
    }
}

impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    ///
    /// compare with other dbg (as answer) and calculate the number of kmers with same copy number
    ///
    pub fn compare(&self, dbg_true: &Dbg<N, E>) -> CompareResult {
        assert_eq!(self.k(), dbg_true.k());
        let p = self.to_kmer_profile();
        let mut r = CompareResult::default();

        for (node, weight) in dbg_true.nodes() {
            let copy_num_true = weight.copy_num();
            let copy_num = match p.get(weight.kmer()) {
                Some(copy_num) => *copy_num,
                None => 0,
            };
            if copy_num_true == copy_num {
                r.n_true += 1;
            } else {
                r.n_error += 1;
            }
        }

        r
    }
    ///
    /// Assuming linear seq
    ///
    pub fn compare_with_seq<S: Seq>(
        &self,
        dbg_true: &Dbg<N, E>,
        seq_true: &S,
    ) -> CompareWithSeqResult<N::Kmer> {
        assert_eq!(self.k(), dbg_true.k());
        let p_self = self.to_kmer_profile();
        let p_true = dbg_true.to_kmer_profile();

        let results: Vec<_> = linear_sequence_to_kmers(seq_true.as_ref(), self.k())
            .map(|kmer| {
                let copy_num_self = match p_self.get(&kmer) {
                    Some(copy_num) => *copy_num,
                    None => 0,
                };
                let copy_num_true = match p_true.get(&kmer) {
                    Some(copy_num) => *copy_num,
                    None => 0,
                };
                CompareWithSeqResultForBase {
                    kmer,
                    copy_num_true,
                    copy_num_self,
                }
            })
            .collect();

        CompareWithSeqResult(results)
    }
    ///
    /// Check k-mer existence (i.e. copy num > 0) of the seq.
    ///
    pub fn check_kmer_existence_with_seqs<T>(&self, seqs: T) -> KmerExistenceResult<N::Kmer>
    where
        T: IntoIterator,
        T::Item: Seq,
    {
        let counts = self.to_kmer_profile();
        let mut r = KmerExistenceResult::new();
        for seq in seqs {
            for kmer in linear_sequence_to_kmers(seq.as_ref(), self.k()) {
                let copy_num = match counts.get(&kmer) {
                    Some(copy_num) => *copy_num,
                    None => 0,
                };
                if copy_num > 0 {
                    r.n_exists += 1;
                } else {
                    r.n_not_exists += 1;
                    r.kmers_not_exists.push(kmer);
                }
            }
        }
        r
    }
    ///
    ///
    pub fn kmer_hists_from_seqs<T>(&self, seqs: T) -> KmerHists
    where
        T: IntoIterator + Clone,
        T::Item: Seq,
    {
        // (1) calculate the true copy numbers in seqs
        let hd: HashDbg<N::Kmer> = HashDbg::from_seqs(self.k(), seqs.clone());
        let copy_nums = hd.to_kmer_profile();
        let max_copy_num = copy_nums.values().max().copied().unwrap();

        // (2) collect counts of dbg
        let counts = self.to_kmer_profile();

        let mut hists = KmerHists::new(max_copy_num);

        // (3) count true kmers with copy_num >= 1.
        for seq in seqs {
            for kmer in linear_sequence_to_kmers(seq.as_ref(), self.k()) {
                let copy_num = *copy_nums.get(&kmer).unwrap();
                let count = match counts.get(&kmer) {
                    Some(count) => *count,
                    None => 0,
                };
                hists.get_hist_mut(copy_num).add(count);
            }
        }

        // (4) count false kmers with copy_num == 0.
        // add x0 kmer (= error and not in true seqs) histogram
        for (kmer, &count) in counts.iter() {
            match copy_nums.get(kmer) {
                Some(&copy_num) => {
                    if copy_num == 0 {
                        hists.get_hist_mut(0).add(count);
                    }
                }
                None => {
                    hists.get_hist_mut(0).add(count);
                }
            }
        }

        hists
    }
}

//
// tests
//
#[cfg(test)]
mod tests {
    use super::*;
    // use crate::common::sequence_to_string;
    use crate::dbg::impls::SimpleDbg;
    use crate::kmer::veckmer::VecKmer;

    #[test]
    fn dbg_compare_simple() {
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_seq(8, b"ATCGGATCGATGC");
        let dbg_true: SimpleDbg<VecKmer> = SimpleDbg::from_seq(8, b"ATCGGCTCGATGC");
        let r = dbg.compare(&dbg_true);
        println!("{:?}", r);
        assert_eq!(r.n_true, 12);
        assert_eq!(r.n_error, 8);

        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_seq(8, b"ATCGGATCGATGC");
        let dbg_true: SimpleDbg<VecKmer> = SimpleDbg::from_seq(8, b"ATCGGATCGATGC");
        let r = dbg.compare(&dbg_true);
        println!("{:?}", r);
        assert_eq!(r.n_true, 20);
        assert_eq!(r.n_error, 0);
    }
    #[test]
    fn dbg_compare_with_seq() {
        let s = b"ATCGGATCGATGC";
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_seq(8, s);
        let dbg_true: SimpleDbg<VecKmer> = SimpleDbg::from_seq(8, s);
        let r = dbg.compare_with_seq(&dbg_true, s);
        println!("{}", r);
    }
    #[test]
    fn dbg_compare_check_kmer_existence() {
        // compare with true seq
        let s = b"ATCGGATCGATGC";
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_seq(8, s);
        let r = dbg.check_kmer_existence_with_seqs(&[s]);
        println!("{:?}", r);
        assert_eq!(r.n_exists, 20);
        assert_eq!(r.n_not_exists, 0);
        assert_eq!(r.kmers_not_exists.len(), 0);

        // compare with false seq with additional A in the last
        let s2 = b"ATCGGATCGATGCA";
        let r = dbg.check_kmer_existence_with_seqs(&[s2]);
        println!("{:?}", r);
        assert_eq!(r.n_exists, 13);
        assert_eq!(r.n_not_exists, 8);
        for kmer in r.kmers_not_exists.iter() {
            println!("{}", kmer);
        }
        assert_eq!(
            r.kmers_not_exists,
            vec![
                VecKmer::from_bases(b"TCGATGCA"),
                VecKmer::from_bases(b"CGATGCAn"),
                VecKmer::from_bases(b"GATGCAnn"),
                VecKmer::from_bases(b"ATGCAnnn"),
                VecKmer::from_bases(b"TGCAnnnn"),
                VecKmer::from_bases(b"GCAnnnnn"),
                VecKmer::from_bases(b"CAnnnnnn"),
                VecKmer::from_bases(b"Annnnnnn"),
            ]
        );
    }
    #[test]
    fn dbg_compare_kmer_hists() {
        // [1] compare with true seq
        let s = b"ATCGGATCGATGC";
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_seq(8, s);
        let kh = dbg.kmer_hists_from_seqs(&[s]);
        assert_eq!(kh.max_copy_num(), 1);
        let c0: Vec<_> = kh.get_hist(0).iter().collect();
        println!("{:?}", c0);
        assert_eq!(c0, vec![]);
        let c1: Vec<_> = kh.get_hist(1).iter().collect();
        println!("{:?}", c1);
        assert_eq!(c1, vec![(1, 20)]);
        println!("{}", kh);
        assert_eq!(kh.to_string(), "x0=;x1=1:20");

        // [2] compare with true seq
        let s1 = b"ATCGGATCGATGC".to_vec();
        let s2 = b"GGCTAGCTGAT".to_vec();
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(8, &[&s1, &s2]);
        // 1
        let kh = dbg.kmer_hists_from_seqs(&[&s1, &s2]);
        println!("{}", kh);
        assert_eq!(kh.to_string(), "x0=;x1=1:38");
        // 2
        let kh = dbg.kmer_hists_from_seqs(&[&s1]);
        println!("{}", kh);
        assert_eq!(kh.to_string(), "x0=1:18;x1=1:20");
    }
}
