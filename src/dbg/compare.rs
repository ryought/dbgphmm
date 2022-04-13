//!
//! Comparison methods of two de Bruijn graphs.
//!
use super::dbg::{Dbg, DbgEdge, DbgNode};

///
/// Result of comparing two Dbgs.
///
#[derive(Clone, Debug, Default)]
pub struct CompareResult {
    pub n_true: usize,
    pub n_error: usize,
}

impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    ///
    ///
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
}
