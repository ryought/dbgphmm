//!
//! Wrapper of VecLike<Prob> (= SparseVec or DenseVec)
//! that can be accessed by NodeIndex
//!
//! TODO Prob can be generalized, so that
//! CopyNums and Freqs can also be the instance of the NodeVec
//!
//! TODO add the edge version of NodeVec
//! that can be accessed by EdgeIndex
//!

use crate::prob::Prob;
use crate::veclike::{DenseVec, SparseVec, VecLike};
pub use petgraph::graph::{EdgeIndex, NodeIndex};
use std::ops::{Index, IndexMut};

#[derive(Debug, Clone)]
pub struct NodeVec<V: VecLike<Prob>>(pub V);
pub type DenseNodeVec = NodeVec<DenseVec<Prob>>;
pub type SparseNodeVec = NodeVec<SparseVec<Prob>>;

impl<V: VecLike<Prob>> Index<NodeIndex> for NodeVec<V> {
    type Output = Prob;
    fn index(&self, index: NodeIndex) -> &Self::Output {
        self.0.get_ref(index.index())
    }
}

impl<V: VecLike<Prob>> IndexMut<NodeIndex> for NodeVec<V> {
    fn index_mut(&mut self, index: NodeIndex) -> &mut Self::Output {
        self.0.get_mut_ref(index.index())
    }
}

impl<V: VecLike<Prob>> NodeVec<V> {
    pub fn new(n_nodes: usize, default_value: Prob) -> Self {
        NodeVec(V::new(n_nodes, default_value))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nodevec() {
        let mut v: SparseNodeVec = SparseNodeVec::new(10, Prob::from_prob(0.0));
        v[NodeIndex::new(0)] = Prob::from_prob(0.1112);
        v[NodeIndex::new(3)] = Prob::from_prob(0.91828);
        println!("{:?}", v);
    }
}
