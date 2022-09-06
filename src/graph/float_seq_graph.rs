//!
//! # `FloatSeqGraph`
//!
//! `SeqGraph` with float-valued copy numbers.
//! It can be created from `FloatDbg`.
//!
//! This file defines conversion from FloatSeqGraph into PHMM(PModel).
//! FloatSeqGraph is petgraph DiGraph whose nodes have a float-valued copy numbers and a base.
//!
use crate::common::NULL_BASE;
use crate::dbg::float::CopyDensity;
use crate::hmmv2::common::{PEdge, PModel, PNode};
use crate::hmmv2::params::PHMMParams;
use crate::prob::Prob;
use derive_new::new;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
use petgraph::Direction;

pub trait FloatSeqNode {
    /// the (float/real-valued) copy number of this node
    fn copy_density(&self) -> CopyDensity;
    /// corresponding base of this node
    fn base(&self) -> u8;
    /// This node is emittable or not?
    /// i.e. self.base != 'X'?
    fn is_emittable(&self) -> bool {
        self.base() != NULL_BASE
    }
    /// (effective) copy_density
    /// (if emittable and otherwise 0.0)
    fn emittable_copy_density(&self) -> CopyDensity {
        if self.is_emittable() {
            self.copy_density()
        } else {
            0.0
        }
    }
}

pub trait FloatSeqEdge {
    /// the (float/real-valued) copy number of this edge if assigned
    fn copy_density(&self) -> Option<CopyDensity>;
}

pub trait FloatSeqGraph {
    /// FloatSeqNode -> PHMM node (PNode)
    ///
    /// PHMMNode has emission and init_prob
    ///
    /// ```text
    /// init_prob = c_l / G (if emittable)
    ///             0       (otherwise)
    /// ```
    fn to_phmm_node(&self, node: NodeIndex, total_copy_density: CopyDensity) -> PNode;
    /// FloatSeqEdge -> PHMM edge (PEdge)
    ///
    /// PHMMEdge has trans_prob
    ///
    /// ```text
    /// trans_prob = (child density) / (sum of parent's childs' density)
    /// ```
    fn to_phmm_edge(&self, edge: EdgeIndex) -> PEdge;
    /// FloatSeqGraph (with FloatSeqNode and FloatSeqEdge) -> PHMM model (PModel)
    fn to_phmm(&self, param: PHMMParams) -> PModel;
    // helper functions
    /// calculate the sum of copy_density of all emittable nodes
    fn total_emittable_copy_density(&self) -> CopyDensity;
    /// calculate the sum of copy_density of all emittable child nodes of the given node
    fn total_emittable_child_copy_density(&self, node: NodeIndex) -> CopyDensity;
}

impl<N: FloatSeqNode, E: FloatSeqEdge> FloatSeqGraph for DiGraph<N, E> {
    fn to_phmm_node(&self, node: NodeIndex, total_copy_density: CopyDensity) -> PNode {
        let node_weight = self.node_weight(node).unwrap();
        let init_prob = if node_weight.is_emittable() {
            Prob::from_prob(node_weight.copy_density()) / Prob::from_prob(total_copy_density)
        } else {
            Prob::from_prob(0.0)
        };
        PNode::new(
            0, // fill the copy_num field a dummy value zero
            init_prob,
            node_weight.is_emittable(),
            node_weight.base(),
        )
    }
    fn to_phmm_edge(&self, edge: EdgeIndex) -> PEdge {
        let (parent, child) = self.edge_endpoints(edge).unwrap();
        let edge_weight = self.edge_weight(edge).unwrap();
        let child_weight = self.node_weight(child).unwrap();
        match edge_weight.copy_density() {
            Some(copy_num) => {
                panic!("FloatSeqEdge with copy density is not supported yet")
            }
            None => {
                // there is no copy num assigned to the edge
                let total_child_copy_density = self.total_emittable_child_copy_density(parent);
                let trans_prob = if child_weight.is_emittable() && total_child_copy_density > 0.0 {
                    Prob::from_prob(child_weight.copy_density())
                        / Prob::from_prob(total_child_copy_density)
                } else {
                    Prob::from_prob(0.0)
                };
                PEdge::new(trans_prob)
            }
        }
    }
    fn to_phmm(&self, param: PHMMParams) -> PModel {
        let total_copy_density = self.total_emittable_copy_density();
        let graph = self.map(
            |v, _| self.to_phmm_node(v, total_copy_density),
            |e, _| self.to_phmm_edge(e),
        );
        PModel { param, graph }
    }
    //
    // helper functions
    //
    fn total_emittable_copy_density(&self) -> CopyDensity {
        self.node_indices()
            .map(|v| self.node_weight(v).unwrap().emittable_copy_density())
            .sum()
    }
    fn total_emittable_child_copy_density(&self, node: NodeIndex) -> CopyDensity {
        self.neighbors_directed(node, Direction::Outgoing)
            .map(|child| self.node_weight(child).unwrap().emittable_copy_density())
            .sum()
    }
}

//
// minimal FloatSeqGraph implementations
// for debugging
//
/// minimal `FloatSeqNode` implementation
#[derive(Debug, Clone, new)]
pub struct SimpleFloatSeqNode {
    copy_density: CopyDensity,
    base: u8,
}
impl FloatSeqNode for SimpleFloatSeqNode {
    fn copy_density(&self) -> CopyDensity {
        self.copy_density
    }
    fn base(&self) -> u8 {
        self.base
    }
}
/// minimal `FloatSeqEdge` implementation
#[derive(Debug, Clone, new)]
pub struct SimpleFloatSeqEdge {}
impl FloatSeqEdge for SimpleFloatSeqEdge {
    fn copy_density(&self) -> Option<CopyDensity> {
        None
    }
}
pub fn mock_simple() -> DiGraph<SimpleFloatSeqNode, SimpleFloatSeqEdge> {
    let mut g = DiGraph::new();
    let va1 = g.add_node(SimpleFloatSeqNode::new(1.0, b'A'));
    let va2 = g.add_node(SimpleFloatSeqNode::new(1.0, b'C'));
    let vb1 = g.add_node(SimpleFloatSeqNode::new(1.0, b'T'));
    let vb2 = g.add_node(SimpleFloatSeqNode::new(1.0, b'G'));
    let vc1 = g.add_node(SimpleFloatSeqNode::new(1.0, b'G'));
    let vc2 = g.add_node(SimpleFloatSeqNode::new(1.0, b'A'));
    g.add_edge(va1, va2, SimpleFloatSeqEdge::new());
    g.add_edge(va2, vb1, SimpleFloatSeqEdge::new());
    g.add_edge(va2, vc1, SimpleFloatSeqEdge::new());
    g.add_edge(vb1, vb2, SimpleFloatSeqEdge::new());
    g.add_edge(vc1, vc2, SimpleFloatSeqEdge::new());
    g
}

//
// tests (with SimpleFloatSeqNode)
//
#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::dot::Dot;

    #[test]
    fn float_seq_graph_simple() {
        let g = mock_simple();
        println!("{:?}", Dot::with_config(&g, &[]));
        let phmm = g.to_phmm(PHMMParams::zero_error());
        println!("{}", phmm);
        // TODO assert the converted phmm has correct init_prob and trans_prob
    }
}
