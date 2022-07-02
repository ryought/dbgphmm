//!
//! Comparison methods of two de Bruijn graphs.
//!
use super::dbg::{Dbg, DbgEdge, DbgNode};
use crate::common::{CopyNum, Seq, SeqStyle, Sequence};
use crate::dbg::hashdbg_v2::HashDbg;
use crate::e2e::Dataset;
use crate::hist::Hist;
use crate::kmer::common::linear_sequence_to_kmers;
use crate::kmer::{KmerLike, NullableKmer};
use crate::prob::Prob;
use itertools::Itertools;

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

impl<K: KmerLike> std::fmt::Display for KmerExistenceResult<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "n_exists={};n_not_exists={}({});",
            self.n_exists,
            self.n_not_exists,
            self.kmers_not_exists.iter().format(","),
        )
    }
}

// kmer classification
#[derive(Clone, Debug, Default)]
pub struct KmerCounts {
    pub n_normal: usize,
    pub n_has_null: usize,
}

impl KmerCounts {
    pub fn add<K: KmerLike>(&mut self, kmer: &K) {
        if kmer.has_null() {
            self.n_has_null += 1;
        } else {
            self.n_normal += 1;
        }
    }
}

impl std::fmt::Display for KmerCounts {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}/{}", self.n_normal, self.n_has_null)
    }
}

///
///
#[derive(Clone, Debug, Default)]
pub struct KmerClassificationResult {
    /// FalsePositive: marked as >0x mistakenly
    pub n_false_kmer_in_dbg: KmerCounts,
    /// TrueNegative: marked as 0x correctly
    pub n_false_kmer_not_in_dbg: KmerCounts,
    /// TruePositive: marked as >0x correctly
    pub n_true_kmer_in_dbg: KmerCounts,
    /// FalseNegative: marked as 0x mistakenly
    pub n_true_kmer_not_in_dbg: KmerCounts,
    ///
    pub n_true_kmer_not_in_reads: KmerCounts,
}

impl std::fmt::Display for KmerClassificationResult {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "TP={};FP={};TN={};FN={};XX={};",
            self.n_true_kmer_in_dbg,
            self.n_false_kmer_in_dbg,
            self.n_false_kmer_not_in_dbg,
            self.n_true_kmer_not_in_dbg,
            self.n_true_kmer_not_in_reads,
        )
    }
}

// kmer histogram

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
    ///
    /// total count of missed kmers,
    /// that is a kmer whose true copy number is 1 or more, but
    /// is not in the dbg. (true-negative)
    ///
    pub fn n_missed_kmers(&self) -> usize {
        (1..=self.max_copy_num())
            .map(|copy_num| self.get_hist(copy_num).get(0))
            .sum()
    }
}

impl std::fmt::Display for KmerHists {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        (0..=self.max_copy_num()).try_for_each(|copy_num| {
            let hist = self.get_hist(copy_num);
            if !hist.is_empty() {
                write!(f, "x{}={};", copy_num, hist)
            } else {
                write!(f, "")
            }
        })
    }
}

//
// CompressionBenchResult
//

///
/// result of compression benchmark
///
pub struct CompressionBenchResult<K: KmerLike> {
    ///
    /// P(R|G)
    ///
    full_prob: Prob,
    ///
    /// genome size of dbg
    ///
    genome_size: CopyNum,
    ///
    ///
    ///
    kmer_existence: KmerExistenceResult<K>,
    ///
    /// histogram of kmer
    ///
    kmer_hists: KmerHists,
    ///
    /// kmer classification result (TP/TN/FP/FN)
    ///
    kmer_classification: KmerClassificationResult,
}

impl<K: KmerLike> std::fmt::Display for CompressionBenchResult<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            self.full_prob.to_log_value(),
            self.genome_size,
            self.kmer_existence,
            self.kmer_classification,
            self.kmer_hists,
        )
    }
}

//
// compare methods for dbg
//

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
    ///
    ///
    pub fn kmer_classification_result<G, R>(&self, genome: G, reads: R) -> KmerClassificationResult
    where
        G: IntoIterator,
        G::Item: Seq,
        R: IntoIterator,
        R::Item: Seq,
    {
        let mut r = KmerClassificationResult::default();

        // counts: dbg
        let counts = self.to_kmer_profile();

        // copy_nums: genome
        let hd_g: HashDbg<N::Kmer> = HashDbg::from_seqs(self.k(), genome);
        let copy_nums = hd_g.to_kmer_profile();

        // read_counts: reads
        let hd_r: HashDbg<N::Kmer> = HashDbg::from_seqs(self.k(), reads);
        let read_counts = hd_r.to_kmer_profile();

        // [1] for all kmers in read
        for (kmer, &read_count) in read_counts.iter() {
            let count = counts.get(&kmer).copied().unwrap_or(0);
            let copy_num = copy_nums.get(&kmer).copied().unwrap_or(0);

            if copy_num > 0 {
                // this is true kmer
                if count > 0 {
                    r.n_true_kmer_in_dbg.add(kmer);
                } else {
                    r.n_true_kmer_not_in_dbg.add(kmer);
                }
            } else {
                // this is false kmer
                if count > 0 {
                    r.n_false_kmer_in_dbg.add(kmer);
                } else {
                    r.n_false_kmer_not_in_dbg.add(kmer);
                }
            }
        }

        // [2] for all kmers in genome but not in read
        for (kmer, &copy_num) in copy_nums.iter() {
            let read_count = read_counts.get(&kmer).copied().unwrap_or(0);
            if read_count == 0 {
                r.n_true_kmer_not_in_reads.add(kmer);
            }
        }

        r
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
        // and create hashmap `copy_nums`
        let hd: HashDbg<N::Kmer> = HashDbg::from_seqs(self.k(), seqs.clone());
        let copy_nums = hd.to_kmer_profile();
        let max_copy_num = copy_nums.values().max().copied().unwrap();

        // (2) collect counts of dbg
        // and create hashmap `counts`
        let counts = self.to_kmer_profile();

        let mut hists = KmerHists::new(max_copy_num);

        // (3) count true kmers that is stored in copy_nums
        for (kmer, &copy_num) in copy_nums.iter() {
            let count = match counts.get(&kmer) {
                Some(count) => *count,
                None => 0,
            };
            hists.get_hist_mut(copy_num).add(count);
        }

        // (4) count false kmers that is not in copy_nums but in counts
        // add x0 kmer (= error and not in true seqs) histogram
        for (kmer, &count) in counts.iter() {
            match copy_nums.get(kmer) {
                Some(_) => {}
                None => {
                    hists.get_hist_mut(0).add(count);
                }
            }
        }

        hists
    }
    ///
    /// benchmark the dbg infered by compression algorithm
    /// using the dataset
    ///
    /// * score
    /// * genome_size
    /// * kmer_existence
    /// * kmer_hists
    ///
    pub fn benchmark_compression(&self, dataset: &Dataset) -> CompressionBenchResult<N::Kmer> {
        let p = self.to_full_prob(dataset.phmm_params.clone(), &dataset.reads);
        let kh = self.kmer_hists_from_seqs(&dataset.genome);
        let ke = self.check_kmer_existence_with_seqs(&dataset.genome);
        let kc = self.kmer_classification_result(&dataset.genome, &dataset.reads);

        CompressionBenchResult {
            full_prob: p,
            genome_size: self.genome_size(),
            kmer_existence: ke,
            kmer_hists: kh,
            kmer_classification: kc,
        }
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
        // 1
        let kh = dbg.kmer_hists_from_seqs(&[s]);
        assert_eq!(kh.max_copy_num(), 1);
        let c0: Vec<_> = kh.get_hist(0).iter().collect();
        println!("{:?}", c0);
        assert_eq!(c0, vec![]);
        let c1: Vec<_> = kh.get_hist(1).iter().collect();
        println!("{:?}", c1);
        assert_eq!(c1, vec![(1, 20)]);
        println!("{}", kh);
        assert_eq!(kh.to_string(), "x1=1:20;");
        assert_eq!(kh.n_missed_kmers(), 0);

        // [2] compare with true seq
        let s1 = b"ATCGGATCGATGC".to_vec();
        let s2 = b"GGCTAGCTGAT".to_vec();
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(8, &[&s1, &s2]);
        // 1
        let kh = dbg.kmer_hists_from_seqs(&[&s1, &s2]);
        println!("{}", kh);
        assert_eq!(kh.to_string(), "x1=1:38;");
        assert_eq!(kh.n_missed_kmers(), 0);
        // 2
        let kh = dbg.kmer_hists_from_seqs(&[&s1]);
        println!("{}", kh);
        assert_eq!(kh.to_string(), "x0=1:18;x1=1:20;");
        assert_eq!(kh.n_missed_kmers(), 0);

        // [3] missing seqs
        let s1 = b"ATCGGATCGATGC".to_vec();
        let s2 = b"GGCTAGCTGAT".to_vec();
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(8, &[&s1]);
        let kh = dbg.kmer_hists_from_seqs(&[&s1, &s2]);
        println!("{}", kh);
        assert_eq!(kh.to_string(), "x1=0:18,1:20;");
        assert_eq!(kh.n_missed_kmers(), 18);

        // [4] duplicating
        let s1 = b"ATCGGATCGATGC".to_vec();
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(8, &[&s1, &s1]);
        let kh = dbg.kmer_hists_from_seqs(&[&s1]);
        println!("{}", kh);
        assert_eq!(kh.to_string(), "x1=2:20;");
    }
}
