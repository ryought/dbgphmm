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
}
