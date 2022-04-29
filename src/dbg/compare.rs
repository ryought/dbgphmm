//!
//! Comparison methods of two de Bruijn graphs.
//!
use super::dbg::{Dbg, DbgEdge, DbgNode};
use crate::common::{CopyNum, Seq, SeqStyle, Sequence};
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
    pub fn check_kmer_existence_with_seq<S: Seq>(&self, seq: &S) -> KmerExistenceResult<N::Kmer> {
        let counts = self.to_kmer_profile();
        let mut r = KmerExistenceResult::new();
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
        r
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
        let r = dbg.check_kmer_existence_with_seq(s);
        println!("{:?}", r);
        assert_eq!(r.n_exists, 20);
        assert_eq!(r.n_not_exists, 0);
        assert_eq!(r.kmers_not_exists.len(), 0);

        // compare with false seq with additional A in the last
        let s2 = b"ATCGGATCGATGCA";
        let r = dbg.check_kmer_existence_with_seq(s2);
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
}
