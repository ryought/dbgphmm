//!
//!
//!
//!
use crate::hmmv2::common::{PHMMEdge, PHMMModel, PHMMNode};
use crate::prob::Prob;
use crate::vector::graph::NodeVec;
use crate::vector::Storage;
use itertools::{chain, Itertools};
use petgraph::graph::NodeIndex;

impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    pub fn to_all_nodes(&self) -> Vec<NodeIndex> {
        self.graph.node_indices().collect()
    }
    pub fn to_childs(&self, nodes: &[NodeIndex]) -> Vec<NodeIndex> {
        nodes
            .iter()
            .flat_map(|&node| self.childs(node).map(|(_, child, _)| child))
            .unique()
            .collect()
    }
    pub fn to_childs_and_us(&self, nodes: &[NodeIndex]) -> Vec<NodeIndex> {
        nodes
            .iter()
            .flat_map(|&node| self.childs(node).map(|(_, parent, _)| parent))
            .chain(nodes.iter().copied())
            .unique()
            .collect()
    }
    pub fn to_parents(&self, nodes: &[NodeIndex]) -> Vec<NodeIndex> {
        nodes
            .iter()
            .flat_map(|&node| self.parents(node).map(|(_, parent, _)| parent))
            .unique()
            .collect()
    }
    pub fn to_parents_and_us(&self, nodes: &[NodeIndex]) -> Vec<NodeIndex> {
        nodes
            .iter()
            .flat_map(|&node| self.parents(node).map(|(_, parent, _)| parent))
            .chain(nodes.iter().copied())
            .unique()
            .collect()
    }
    pub fn merge(&self, a: &[NodeIndex], b: &[NodeIndex]) -> Vec<NodeIndex> {
        chain!(a.iter().copied(), b.iter().copied())
            .unique()
            .collect()
    }
}
