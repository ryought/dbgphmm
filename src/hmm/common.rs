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

pub trait PHMMNode {
    fn emission(&self) -> u8;
    fn is_emittable(&self) -> bool {
        self.emission() != b'N'
    }
    // not necessary
    // fn copy_num(&self) -> CopyNum;
    fn init_prob(&self) -> Prob;
}

pub trait PHMMEdge {
    // TODO not necesssary
    // fn copy_num(&self) -> CopyNum;
    fn trans_prob(&self) -> Prob;
}

pub struct PHMMGraph<N: PHMMNode, E: PHMMEdge> {
    param: PHMMParams,
    grpah: DiGraph<N, E>,
}

pub type PGraph = PHMMGraph<PNode, PEdge>;

pub struct PNode {
    ///
    copy_num: CopyNum,
    /// initial transition probability
    ///  = (copy_num) / (sum of copy_nums of all nodes)
    init_prob: Prob,
    ///
    is_emittable: bool,
    ///
    emission: u8,
}

impl PHMMNode for PNode {
    fn emission(&self) -> u8 {
        self.emission
    }
    fn is_emittable(&self) -> bool {
        self.is_emittable
    }
    fn init_prob(&self) -> Prob {
        self.init_prob
    }
}

pub struct PEdge {
    ///
    trans_prob: Prob,
    ///
    copy_num: Option<CopyNum>,
}

impl PEdge {
    pub fn new(trans_prob: Prob) -> PEdge {
        PEdge {
            trans_prob,
            copy_num: None,
        }
    }
}

impl PHMMEdge for PEdge {
    fn trans_prob(&self) -> Prob {
        self.trans_prob
    }
}

//
// Forward related
//
pub struct PHMMForwardResult<V: VecLike<Prob>>(Vec<PHMMTable<V>>);

pub const MAX_DEL: usize = 4;

#[derive(Debug, Clone)]
pub struct PHMMTable<V: VecLike<Prob>> {
    /// Match node probability
    m: V,
    /// Ins node probability
    i: V,
    /// Del node probability
    d: V,
    /// Match node in begin state
    mb: Prob,
    /// Ins node in begin state
    ib: Prob,
    /// end state probability
    e: Prob,
}

impl<N: PHMMNode, E: PHMMEdge> PHMMGraph<N, E> {
    ///
    fn forward<V: VecLike<Prob>>(&self) -> PHMMForwardResult<V> {
        unimplemented!();
    }
}

// impl PHMMForward for PHMM {}
