//!
//! Definition of Node-centric PHMM
//!
use super::params::PHMMParams;
use crate::common::CopyNum;
use crate::prob::Prob;
use crate::veclike::{DenseVec, VecLike};
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

//
// traits
//

trait PHMMNode {
    fn emission(&self) -> u8;
    fn is_emittable(&self) -> bool {
        self.emission() != b'N'
    }
    fn copy_num(&self) -> CopyNum;
    fn init_prob(&self) -> Prob;
}

trait PHMMEdge {
    fn copy_num(&self) -> CopyNum;
}

struct PHMM<N: PHMMNode, E: PHMMEdge> {
    param: PHMMParams,
    grpah: DiGraph<N, E>,
}

struct PNode {
    ///
    copy_num: CopyNum,
    /// initial transition probability
    ///  = (copy_num) / (sum of copy_nums of all nodes)
    init_prob: Prob,
    ///
    is_emitable: bool,
}

struct PEdge {
    ///
    trans_prob: Prob,
    ///
    copy_num: Option<CopyNum>,
}

// impl PHMM {}

//
// Forward related
//
struct PHMMForwardResult<V: VecLike<Prob>>(Vec<PHMMTable<V>>);

pub const MAX_DEL: usize = 4;

#[derive(Debug, Clone)]
pub struct PHMMTable<V: VecLike<Prob>> {
    m: V,
    i: V,
    d: V,
    mb: Prob,
    ib: Prob,
    e: Prob,
}

trait PHMMForward {
    fn forward<V: VecLike<Prob>>(&self) -> PHMMForwardResult<V> {
        unimplemented!();
    }
}

// impl PHMMForward for PHMM {}
