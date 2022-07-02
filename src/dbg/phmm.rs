//!
//! dbg as a seqgraph and phmm
//!
use super::dbg::{Dbg, DbgEdge, DbgNode};
use crate::common::{CopyNum, Seq};
use crate::graph::seq_graph::{SeqEdge, SeqGraph, SeqNode};
use crate::hmmv2::common::PModel;
use crate::hmmv2::params::PHMMParams;
use crate::prob::Prob;
use rayon::prelude::*;

impl<N: DbgNode> SeqNode for N {
    fn copy_num(&self) -> CopyNum {
        self.copy_num()
    }
    fn base(&self) -> u8 {
        self.emission()
    }
}

impl<E: DbgEdge> SeqEdge for E {
    fn copy_num(&self) -> Option<CopyNum> {
        self.copy_num()
    }
}

impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    /// Convert dbg into phmm
    pub fn to_phmm(&self, param: PHMMParams) -> PModel {
        self.graph.to_phmm(param)
    }
    /// Convert dbg into phmm and calculate full probability
    pub fn to_full_prob<T>(&self, param: PHMMParams, seqs: T) -> Prob
    where
        T: IntoParallelIterator,
        T::Item: Seq,
    {
        let phmm = self.to_phmm(param);
        phmm.to_full_prob_parallel(seqs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dbg::mocks::*;
    #[test]
    fn dbg_as_phmm_simple() {
        let dbg = mock_simple();
        let c = dbg.graph.total_emittable_copy_num();
        println!("c={}", c);
        let phmm = dbg.to_phmm(PHMMParams::default());
        println!("{}", phmm);
    }
}
