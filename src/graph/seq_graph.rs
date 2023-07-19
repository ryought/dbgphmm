//!
//! `SeqGraph` is `DiGraph<N: SeqNode, E: SeqEdge>`
//! seq with copy numbers
//!

use crate::common::{CopyNum, NULL_BASE};
use crate::graph::genome_graph::GenomeGraphPos;
use crate::hmmv2::common::{PEdge, PModel, PNode};
use crate::hmmv2::params::PHMMParams;
use crate::prob::Prob;
use petgraph::dot::Dot;
use petgraph::graph::DiGraph;
pub use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::visit::{IntoNodeReferences, NodeRef};
pub use petgraph::Direction;

pub trait SeqGraph {
    /// calculate the sum of copy numbers
    /// of all emittable nodes
    fn total_emittable_copy_num(&self, min_copy_num: CopyNum) -> CopyNum;
    /// calculate the sum of copy numbers
    /// of all emittable childs of the given node
    fn total_emittable_child_copy_nums(&self, node: NodeIndex, min_copy_num: CopyNum) -> CopyNum;
    ///
    /// SeqGraph has consistent copy numbers on nodes?
    ///
    /// * for all nodes, sum of copy numbers of in-edges and out-edges are the same.
    ///
    fn node_copy_nums_is_consistent(&self) -> bool;
    ///
    /// SeqGraph has consistent copy numbers on edges?
    ///
    /// * for all nodes, all of out-edges are either with-copy-numbers or without-copy-numbers
    /// * if copy numbers are set on edges, the sum of copy numbers on edges should be equal to the
    /// copy number of the node.
    ///
    fn edge_copy_nums_is_consistent(&self) -> bool;
    ///
    /// Check if all edges out-going from this node
    /// has its own copy_nums
    ///
    fn all_out_edges_has_copy_nums(&self, node: NodeIndex) -> bool;
    ///
    /// Check if no edges out-going from this node
    /// has its own copy nums
    ///
    fn no_out_edges_has_copy_nums(&self, node: NodeIndex) -> bool;
    ///
    /// Sum of out-edge copy_nums
    ///
    fn sum_out_edge_copy_nums(&self, node: NodeIndex) -> CopyNum;
    //
    // normal phmm
    //
    /// Convert Node in SimpleSeqGraph into phmm node
    ///
    fn to_phmm_node(
        &self,
        node: NodeIndex,
        total_copy_num: CopyNum,
        min_copy_num: CopyNum,
    ) -> PNode;
    /// Convert Edge in SimpleSeqGraph into phmm edge
    ///
    /// For each node,
    /// * if all outgoing edge do not has copynum,
    ///
    /// * if all outgoing edge have copynum,
    ///
    fn to_phmm_edge(&self, edge: EdgeIndex, min_copy_num: CopyNum) -> PEdge;
    /// convert SimpleSeqGraph to PHMM by ignoreing the edge copy numbers
    ///
    fn to_phmm(&self, param: PHMMParams) -> PModel;
    //
    // uniform phmm
    //
    /// Create node for uniform PHMM
    ///
    /// initial probability of node `v` =
    ///   `1 / n_emittable_nodes` (if `v` is emittable)
    ///   `0` (otherwise)
    ///
    fn to_uniform_phmm_node(&self, node: NodeIndex, n_emittable_nodes: usize) -> PNode;
    /// Create edge for uniform PHMM
    ///
    /// transition probability of edge `(s,t)` =
    ///   `1 / n_emittable_childs` (if `t` is emittable)
    ///   `0` (otherwise)
    ///
    fn to_uniform_phmm_edge(&self, edge: EdgeIndex) -> PEdge;
    /// Craete uniform PHMM for creating hint of reads
    ///
    ///
    fn to_uniform_phmm(&self, param: PHMMParams) -> PModel;
    //
    // non-zero phmm
    //
    ///
    /// Non-zero PHMM
    ///
    /// Set init_prob and trans_prob according to node copy numbers
    /// but copy_num is clamped to 1 if the node is emittable.
    /// Use to create mapping.
    ///
    fn to_non_zero_phmm(&self, param: PHMMParams) -> PModel;
}

impl<N: SeqNode, E: SeqEdge> SeqGraph for DiGraph<N, E> {
    /// calculate the sum of copy numbers
    /// of all emittable nodes
    fn total_emittable_copy_num(&self, min_copy_num: CopyNum) -> CopyNum {
        self.node_indices()
            .map(|v| {
                let vw = self.node_weight(v).unwrap();
                if vw.is_emittable() {
                    vw.copy_num().max(min_copy_num)
                } else {
                    0
                }
            })
            .sum()
    }
    /// calculate the sum of copy numbers
    /// of all emittable childs of the given node
    fn total_emittable_child_copy_nums(&self, node: NodeIndex, min_copy_num: CopyNum) -> CopyNum {
        self.neighbors_directed(node, Direction::Outgoing)
            .map(|child| {
                let child_weight = self.node_weight(child).unwrap();
                if child_weight.is_emittable() {
                    child_weight.copy_num().max(min_copy_num)
                } else {
                    0
                }
            })
            .sum()
    }
    fn node_copy_nums_is_consistent(&self) -> bool {
        // TODO
        true
    }
    fn edge_copy_nums_is_consistent(&self) -> bool {
        self.node_indices().all(|v| {
            let vw = self.node_weight(v).unwrap();
            self.edges_directed(v, Direction::Outgoing)
                .all(|e| e.weight().copy_num().is_some())
        })
    }
    fn all_out_edges_has_copy_nums(&self, node: NodeIndex) -> bool {
        self.edges_directed(node, Direction::Outgoing)
            .all(|e| e.weight().copy_num().is_some())
    }
    fn no_out_edges_has_copy_nums(&self, node: NodeIndex) -> bool {
        self.edges_directed(node, Direction::Outgoing)
            .all(|e| e.weight().copy_num().is_none())
    }
    fn sum_out_edge_copy_nums(&self, node: NodeIndex) -> CopyNum {
        self.edges_directed(node, Direction::Outgoing)
            .map(|e| e.weight().copy_num().unwrap())
            .sum()
    }
    /// Convert Node in SimpleSeqGraph into phmm node
    fn to_phmm_node(
        &self,
        node: NodeIndex,
        total_copy_num: CopyNum,
        min_copy_num: CopyNum,
    ) -> PNode {
        let node_weight = self.node_weight(node).unwrap();
        let init_prob = if node_weight.is_emittable() {
            let copy_num = node_weight.copy_num().max(min_copy_num);
            Prob::from_prob(copy_num as f64) / Prob::from_prob(total_copy_num as f64)
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
    fn to_phmm_edge(&self, edge: EdgeIndex, min_copy_num: CopyNum) -> PEdge {
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
                let total_child_copy_num =
                    self.total_emittable_child_copy_nums(parent, min_copy_num);
                let trans_prob = if child_weight.is_emittable() && total_child_copy_num > 0 {
                    let copy_num = child_weight.copy_num().max(min_copy_num);
                    Prob::from_prob(copy_num as f64 / total_child_copy_num as f64)
                } else {
                    Prob::from_prob(0.0)
                };
                PEdge::new(trans_prob)
            }
        }
    }
    /// convert SimpleSeqGraph to PHMM by ignoreing the edge copy numbers
    fn to_phmm(&self, param: PHMMParams) -> PModel {
        let min_copy_num = 0;
        let total_copy_num = self.total_emittable_copy_num(min_copy_num);
        let graph = self.map(
            // node converter
            |v, _| self.to_phmm_node(v, total_copy_num, min_copy_num),
            // edge converter
            |e, _| self.to_phmm_edge(e, min_copy_num),
        );
        PModel { param, graph }
    }
    fn to_uniform_phmm_node(&self, node: NodeIndex, n_emittable_nodes: usize) -> PNode {
        let node_weight = self.node_weight(node).unwrap();
        let init_prob = if node_weight.is_emittable() {
            Prob::one() / Prob::from_prob(n_emittable_nodes as f64)
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
    fn to_uniform_phmm_edge(&self, edge: EdgeIndex) -> PEdge {
        let (parent, child) = self.edge_endpoints(edge).unwrap();
        let n_emittable_childs = self
            .neighbors_directed(parent, Direction::Outgoing)
            .filter(|&v| self.node_weight(v).unwrap().is_emittable())
            .count();
        let is_emittable = self.node_weight(child).unwrap().is_emittable();
        let trans_prob = if is_emittable {
            Prob::one() / Prob::from_prob(n_emittable_childs as f64)
        } else {
            Prob::from_prob(0.0)
        };
        PEdge::new(trans_prob)
    }
    fn to_uniform_phmm(&self, param: PHMMParams) -> PModel {
        let n_emittable_nodes = self
            .node_references()
            .filter(|&v| v.weight().is_emittable())
            .count();
        let graph = self.map(
            |v, _| self.to_uniform_phmm_node(v, n_emittable_nodes),
            |e, _| self.to_uniform_phmm_edge(e),
        );
        PModel { param, graph }
    }
    fn to_non_zero_phmm(&self, param: PHMMParams) -> PModel {
        let min_copy_num = 1;
        let total_copy_num = self.total_emittable_copy_num(min_copy_num);
        let graph = self.map(
            // node converter
            |v, _| self.to_phmm_node(v, total_copy_num, min_copy_num),
            // edge converter
            |e, _| self.to_phmm_edge(e, min_copy_num),
        );
        PModel { param, graph }
    }
}

/// a trait that should be satisfyed by nodes in SeqGraph
pub trait SeqNode {
    ///
    /// the copy number of this node
    ///
    fn copy_num(&self) -> CopyNum;
    ///
    /// corresponding base of this node
    ///
    fn base(&self) -> u8;
    ///
    /// This node is emittable or not?
    /// i.e. self.base != 'X'?
    ///
    fn is_emittable(&self) -> bool {
        self.base() != NULL_BASE
    }
}

/// a trait that should be satisfyed by edges in SeqGraph
pub trait SeqEdge {
    ///
    /// the copy number of this edge if assigned
    ///
    fn copy_num(&self) -> Option<CopyNum>;
}

//
// minimum implementations
//

///
/// SimpleSeqGraph is a sequence graph whose node has its own copy number.
///
pub type SimpleSeqGraph = DiGraph<SimpleSeqNode, SimpleSeqEdge>;

/// Get a vector of start point node index
pub fn get_start_points(g: &DiGraph<SimpleSeqNode, SimpleSeqEdge>) -> Vec<NodeIndex> {
    g.node_indices()
        .filter(|&v| g.node_weight(v).unwrap().is_start_point())
        .collect()
}

/// Find a node in seqgraph
/// with given GenomeGraphPos and revcomp info
pub fn find_node_from_source_pos(
    g: &DiGraph<SimpleSeqNode, SimpleSeqEdge>,
    source: GenomeGraphPos,
    is_revcomp: bool,
) -> Option<NodeIndex> {
    g.node_indices().find(|&v| {
        let weight = g.node_weight(v).unwrap();
        weight.source == source && weight.is_revcomp == is_revcomp
    })
}

pub struct SimpleSeqNode {
    copy_num: CopyNum,
    base: u8,
    is_start_point: bool,
    is_revcomp: bool,
    source: GenomeGraphPos,
}

impl SimpleSeqNode {
    pub fn new(
        copy_num: CopyNum,
        base: u8,
        is_start_point: bool,
        is_revcomp: bool,
        source: GenomeGraphPos,
    ) -> SimpleSeqNode {
        SimpleSeqNode {
            copy_num,
            base,
            is_start_point,
            is_revcomp,
            source,
        }
    }
    pub fn is_start_point(&self) -> bool {
        self.is_start_point
    }
    pub fn source(&self) -> GenomeGraphPos {
        self.source
    }
    pub fn is_revcomp(&self) -> bool {
        self.is_revcomp
    }
}

impl SeqNode for SimpleSeqNode {
    fn copy_num(&self) -> CopyNum {
        self.copy_num
    }
    fn base(&self) -> u8 {
        self.base
    }
}

pub struct SimpleSeqEdge {
    copy_num: Option<CopyNum>,
}

impl SimpleSeqEdge {
    pub fn new(copy_num: Option<CopyNum>) -> SimpleSeqEdge {
        SimpleSeqEdge { copy_num }
    }
}

impl SeqEdge for SimpleSeqEdge {
    fn copy_num(&self) -> Option<CopyNum> {
        self.copy_num
    }
}

//
// Display
//

impl std::fmt::Display for SimpleSeqNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{} (x{}) (s={}, revcomp={}, source={})",
            self.base() as char,
            self.copy_num(),
            self.is_start_point(),
            self.is_revcomp,
            self.source,
        )
    }
}

impl std::fmt::Display for SimpleSeqEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self.copy_num() {
            Some(copy_num) => {
                write!(f, "x{}", copy_num)
            }
            None => {
                write!(f, "")
            }
        }
    }
}

//
// mock constructors
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{ei, ni};
    use crate::graph::mocks::{mock_crossing, mock_linear};
    use crate::hmmv2::params::PHMMParams;
    use crate::prob::p;
    #[test]
    fn trait_test() {
        let sg = mock_linear().to_seq_graph();
    }
    #[test]
    fn seq_graph_edge_copy_num() {
        let m1 = mock_crossing(false);
        let g1 = m1.to_seq_graph().to_phmm(PHMMParams::default());
        println!("{}", m1);
        let m2 = mock_crossing(true);
        let g2 = m2.to_seq_graph().to_phmm(PHMMParams::default());
        println!("{}", m2);

        let ra = b"ATTAGGAGCA";
        let rb = b"ATTAGGAGCAGCTGATAGGG";

        let o1a = g1.run(ra);
        let o1b = g1.run(rb);
        println!("{}", o1a.to_full_prob_forward());
        println!("{}", o1b.to_full_prob_backward());
        // assert similarity between forward and backward
        assert_abs_diff_eq!(
            o1a.to_full_prob_forward(),
            o1a.to_full_prob_backward(),
            epsilon = 0.1
        );
        assert_abs_diff_eq!(
            o1b.to_full_prob_forward(),
            o1b.to_full_prob_backward(),
            epsilon = 0.1
        );

        let o2a = g2.run(ra);
        let o2b = g2.run(rb);
        println!("{}", o2a.to_full_prob_forward());
        println!("{}", o2b.to_full_prob_backward());
        // assert similarity between forward and backward
        assert_abs_diff_eq!(
            o2a.to_full_prob_forward(),
            o2a.to_full_prob_backward(),
            epsilon = 0.1
        );
        assert_abs_diff_eq!(
            o2b.to_full_prob_forward(),
            o2b.to_full_prob_backward(),
            epsilon = 0.1
        );

        // assert crossing
        assert_abs_diff_eq!(
            o1a.to_full_prob_forward(),
            o2a.to_full_prob_forward(),
            epsilon = 0.1
        );
        // assert crossing
        assert!(o1b.to_full_prob_forward().to_log_value() > -17.0);
        assert!(o2b.to_full_prob_forward().to_log_value() < -39.0);

        println!("{}", g1);
        let efo1b = o1b.to_edge_freqs(&g1, rb);
        println!("{} {}", efo1b.len(), efo1b);
        assert!(efo1b[ei(36)] < 0.0001);
        assert!(efo1b[ei(37)] > 0.9);
        assert!(efo1b[ei(38)] < 0.0001);
        assert!(efo1b[ei(39)] < 0.0001);
        let efo2b = o2b.to_edge_freqs(&g2, rb);
        println!("{} {}", efo2b.len(), efo2b);
        assert!(efo2b[ei(37)] == 0.0);
        assert!(efo2b[ei(38)] == 0.0);
    }
}
