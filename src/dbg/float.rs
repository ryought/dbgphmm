//!
//! de bruijn graph with float (real-valued) copy numbers
//!
use super::dbg::{Dbg, DbgEdge, DbgNode};
use crate::common::CopyNum;
use crate::graph::float_seq_graph::{FloatSeqEdge, FloatSeqNode};
use crate::kmer::kmer::{Kmer, KmerLike};

/// `CopyDensity` = f64
/// Float valued copy number
pub type CopyDensity = f64;

/// de bruijn graph with float (real-valued) copy numbers
pub type FloatDbg<K> = Dbg<FloatDbgNode<K>, FloatDbgEdge>;

/// node struct of FloatDbg
#[derive(Clone)]
pub struct FloatDbgNode<K: KmerLike> {
    kmer: K,
    copy_density: CopyDensity,
}

impl<K: KmerLike> DbgNode for FloatDbgNode<K> {
    type Kmer = K;
    fn new(kmer: K, copy_num: CopyNum) -> Self {
        // cast to f64
        let copy_density = copy_num as CopyDensity;
        FloatDbgNode { kmer, copy_density }
    }
    fn kmer(&self) -> &K {
        &self.kmer
    }
    fn copy_num(&self) -> CopyNum {
        self.copy_density.round() as CopyNum
    }
    fn set_copy_num(&mut self, copy_num: CopyNum) {
        // cast to f64 and store it as copy_density
        self.copy_density = copy_num as CopyDensity;
    }
}

/// edge struct of FloatDbg
#[derive(Clone)]
pub struct FloatDbgEdge {
    copy_density: Option<CopyDensity>,
}

impl DbgEdge for FloatDbgEdge {
    fn new(copy_num: Option<CopyNum>) -> Self {
        // cast to f64
        let copy_density = match copy_num {
            Some(copy_num) => Some(copy_num as CopyDensity),
            None => None,
        };
        FloatDbgEdge { copy_density }
    }
    fn copy_num(&self) -> Option<CopyNum> {
        match self.copy_density {
            Some(copy_density) => Some(copy_density.round() as CopyNum),
            None => None,
        }
    }
    fn set_copy_num(&mut self, copy_num: Option<CopyNum>) {
        self.copy_density = match copy_num {
            Some(copy_num) => Some(copy_num as CopyDensity),
            None => None,
        };
    }
}

//
// std::fmt::Display
//
impl<K: KmerLike> std::fmt::Display for FloatDbgNode<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{} (x{})", self.kmer(), self.copy_density)
    }
}
impl std::fmt::Display for FloatDbgEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self.copy_density {
            Some(copy_density) => write!(f, "x{}", copy_density),
            None => write!(f, "x?"),
        }
    }
}

//
// FloatSeqGraph
//
impl<K: KmerLike> FloatSeqNode for FloatDbgNode<K> {
    fn copy_density(&self) -> CopyDensity {
        self.copy_density
    }
    fn base(&self) -> u8 {
        self.emission()
    }
}
impl FloatSeqEdge for FloatDbgEdge {
    fn copy_density(&self) -> Option<CopyDensity> {
        self.copy_density
    }
}
