//!
//! # `FloatSeqGraph`
//!
//! `SeqGraph` with float-valued copy numbers.
//! It can be created from `FloatDbg`.
//!
//! This file defines conversion from FloatSeqGraph into PHMM(PModel).
//! FloatSeqGraph is petgraph DiGraph whose nodes have a float-valued copy numbers and a base.
//!
use crate::common::{CopyNum, NULL_BASE};
use crate::dbg::float::CopyDensity;
use crate::hmmv2::common::{PEdge, PModel, PNode};
use crate::hmmv2::params::PHMMParams;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

pub trait FloatSeqNode {
    /// the (float/real-valued) copy number of this node
    fn copy_density(&self) -> CopyDensity;
    /// corresponding base of this node
    fn base(&self) -> u8;
    // /// This node is emittable or not?
    // /// i.e. self.base != 'X'?
    // fn is_emittable(&self) -> bool {
    //     self.base() != NULL_BASE
    // }
}

pub trait FloatSeqEdge {
    /// the (float/real-valued) copy number of this edge if assigned
    fn copy_density(&self) -> Option<CopyDensity>;
}

pub trait FloatSeqGraph {
    /*
    fn to_phmm_node(&self, node: NodeIndex, total_copy_num: CopyNum) -> PNode;
    fn to_phmm_edge(&self, edge: EdgeIndex) -> PEdge;
    fn to_phmm(&self, param: PHMMParams) -> PModel;
    */
}

impl<N: FloatSeqNode, E: FloatSeqEdge> FloatSeqGraph for DiGraph<N, E> {
    /*
    /// Convert Node in SimpleSeqGraph into phmm node
    fn to_phmm_node(&self, node: NodeIndex, total_copy_num: CopyNum) -> PNode {
        let node_weight = self.node_weight(node).unwrap();
        let init_prob = if node_weight.is_emittable() {
            Prob::from_prob(node_weight.copy_num() as f64) / Prob::from_prob(total_copy_num as f64)
        } else {
            Prob::from_prob(0.0)
        };
        PNode::new(
            node_weight.copy_num(),
            init_prob,
            node_weight.is_emittable(),
            node_weight.base(),
        )
    }
    fn to_phmm_edge(&self, edge: EdgeIndex) -> PEdge {
        let (parent, child) = self.edge_endpoints(edge).unwrap();
        let edge_weight = self.edge_weight(edge).unwrap();
        let child_weight = self.node_weight(child).unwrap();
        match edge_weight.copy_num() {
            Some(copy_num) => {
                // copy num is assigned
                // XXX assume their consistency on the graph
                let parent_weight = self.node_weight(parent).unwrap();
                let parent_copy_num = parent_weight.copy_num();
                let trans_prob = if child_weight.is_emittable() && copy_num > 0 {
                    assert!(parent_copy_num > 0);
                    Prob::from_prob(copy_num as f64 / parent_copy_num as f64)
                } else {
                    Prob::from_prob(0.0)
                };
                PEdge::new(trans_prob)
            }
            None => {
                // there is no copy num assigned to the edge
                let total_child_copy_num = self.total_emittable_child_copy_nums(parent);
                let trans_prob = if child_weight.is_emittable() && total_child_copy_num > 0 {
                    Prob::from_prob(child_weight.copy_num() as f64 / total_child_copy_num as f64)
                } else {
                    Prob::from_prob(0.0)
                };
                PEdge::new(trans_prob)
            }
        }
    }
    fn to_phmm(&self, param: PHMMParams) -> PModel {
        let total_copy_num = self.total_emittable_copy_num();
        let graph = self.map(
            // node converter
            |v, _| self.to_phmm_node(v, total_copy_num),
            // edge converter
            |e, _| self.to_phmm_edge(e),
        );
        PModel { param, graph }
    }
    */
}
