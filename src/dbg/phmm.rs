//!
//! dbg as a seqgraph and phmm
//!
use super::dbg::{Dbg, DbgEdge, DbgNode};
use crate::common::CopyNum;
use crate::graph::seq_graph::{SeqEdge, SeqGraph, SeqNode};

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
