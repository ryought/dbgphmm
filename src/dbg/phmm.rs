//!
//! dbg as a seqgraph and phmm
//!
use super::dbg::{Dbg, DbgEdge, DbgNode};
use crate::common::CopyNum;
use crate::graph::seq_graph::{SeqEdge, SeqGraph, SeqNode};
use crate::hmm::params::PHMMParams;
use crate::hmmv2::common::PModel;

impl<N> SeqNode for N
where
    N: DbgNode,
{
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
    fn to_phmm(&self, param: PHMMParams) -> PModel {
        self.graph.to_phmm(param)
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
