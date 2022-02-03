//!
//! Definition of Node-centric PHMM
//!
use super::params::PHMMParams;
use crate::common::CopyNum;
use crate::prob::Prob;
use crate::veclike::{DenseVec, VecLike};
use petgraph::dot::Dot;
use petgraph::graph::DiGraph;
pub use petgraph::graph::{EdgeIndex, NodeIndex};

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

pub struct PHMMModel<N: PHMMNode, E: PHMMEdge> {
    pub param: PHMMParams,
    pub graph: DiGraph<N, E>,
}

pub type PGraph = DiGraph<PNode, PEdge>;
pub type PModel = PHMMModel<PNode, PEdge>;

impl<N: PHMMNode + std::fmt::Display, E: PHMMEdge + std::fmt::Display> PHMMModel<N, E> {
    pub fn to_dot(&self) {
        println!("{}", Dot::with_config(&self.graph, &[]));
    }
}

//
//
// PNode
//
//

#[derive(Debug, Copy, Clone)]
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

impl PNode {
    pub fn new(copy_num: CopyNum, init_prob: Prob, is_emittable: bool, emission: u8) -> PNode {
        PNode {
            copy_num,
            init_prob,
            is_emittable,
            emission,
        }
    }
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

impl std::fmt::Display for PNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.is_emittable {
            write!(f, "{} (p={})", self.emission as char, self.init_prob)
        } else {
            write!(f, "not_emittable (p={})", self.init_prob)
        }
    }
}

//
//
// PEdge
//
//

#[derive(Debug, Copy, Clone)]
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

impl std::fmt::Display for PEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "p={}", self.trans_prob)
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

impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    fn forward<V: VecLike<Prob>>(&self) -> PHMMForwardResult<V> {
        unimplemented!();
    }
}

// impl PHMMForward for PHMM {}
