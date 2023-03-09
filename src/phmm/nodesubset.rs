//!
//! Subset of nodes in PHMM
//!
//! * ActiveNodes
//! * Hint
//!
use super::PHMM;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

///
/// Subset of nodes
///
#[derive(Clone, Debug)]
pub struct NodeSubset(Vec<NodeIndex>);

pub struct NodeSubsetIterator<'a> {
    index: usize,
    subset: &'a NodeSubset,
}

impl NodeSubset {
    ///
    /// Iterate over nodes
    ///
    pub fn iter<'a>(&'a self) -> NodeSubsetIterator<'a> {
        NodeSubsetIterator {
            index: 0,
            subset: self,
        }
    }
}

impl<'a> Iterator for NodeSubsetIterator<'a> {
    type Item = NodeIndex;
    fn next(&mut self) -> Option<Self::Item> {
        let index = self.index;
        if index < self.subset.0.len() {
            self.index += 1;
            Some(self.subset.0[index])
        } else {
            None
        }
    }
}

impl PHMM {
    ///
    /// Get a NodeSubset of childs of nodes in current NodeSubset
    ///
    pub fn to_child_subsets(&self, nodes: &NodeSubset) -> NodeSubset {
        unimplemented!();
    }
}
