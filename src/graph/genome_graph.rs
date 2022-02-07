//!
//! `GenomeGraph`
//! Node: linear sequence
//! Edge: its adjacency
//!
use super::seq_graph::{SeqGraph, SimpleSeqEdge, SimpleSeqNode};
use crate::common::{CopyNum, Sequence};
use itertools::Itertools;
use petgraph::dot::Dot;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use petgraph::visit::IntoNodeReferences;
use std::collections::HashMap;

/// GenomeGraph
pub struct GenomeGraph(pub DiGraph<GenomeNode, GenomeEdge>);

/// Node in GenomeGraph is a copy-number-assigned sequence fragment
pub struct GenomeNode {
    /// sequence fragment of this node
    seq: Sequence,
    /// copy number of the sequence
    copy_num: CopyNum,
}

impl GenomeNode {
    pub fn new(seq: &[u8], copy_num: CopyNum) -> Self {
        GenomeNode {
            seq: seq.to_vec(),
            copy_num,
        }
    }
}

pub struct GenomeEdge {
    /// copy number of the transition
    /// it assumes random transition if `None`
    copy_num: Option<CopyNum>,
}

impl GenomeEdge {
    pub fn new(copy_num: Option<CopyNum>) -> Self {
        GenomeEdge { copy_num }
    }
}

impl std::fmt::Display for GenomeNode {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let seq = std::str::from_utf8(&self.seq).unwrap();
        write!(f, "{} (x{})", seq, self.copy_num)
    }
}

impl std::fmt::Display for GenomeEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self.copy_num {
            Some(copy_num) => {
                write!(f, "x{}", copy_num)
            }
            None => {
                write!(f, "")
            }
        }
    }
}

impl std::fmt::Display for GenomeGraph {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", Dot::with_config(&self.0, &[]))
    }
}

impl GenomeGraph {
    /// Get the number of nodes in the GenomeGraph
    pub fn node_count(&self) -> usize {
        self.0.node_count()
    }
    /// Get the number of edges in the GenomeGraph
    pub fn edge_count(&self) -> usize {
        self.0.edge_count()
    }
    ///
    /// Convert `GenomeGraph` into `SeqGraph`
    /// split a node containing bases into multiple nodes corresponding each bases
    ///
    pub fn to_seq_graph(&self) -> SeqGraph<SimpleSeqNode, SimpleSeqEdge> {
        let mut graph = DiGraph::new();

        // for each node
        let mut m: HashMap<NodeIndex, (NodeIndex, NodeIndex)> = HashMap::new();
        for (node, node_weight) in self.0.node_references() {
            let seq = &node_weight.seq;
            let copy_num = node_weight.copy_num;

            // add a new node for each bases
            let nodes: Vec<NodeIndex> = seq
                .iter()
                .map(|&base| graph.add_node(SimpleSeqNode::new(copy_num, base)))
                .collect();

            // add edges between adjacent two bases
            let edges: Vec<EdgeIndex> = nodes
                .iter()
                .tuple_windows()
                .map(|(&v0, &v1)| graph.add_edge(v0, v1, SimpleSeqEdge::new(None)))
                .collect();

            // store the first/last node
            m.insert(node, (*nodes.first().unwrap(), *nodes.last().unwrap()));
        }

        // for each node
        for er in self.0.edge_references() {
            let source = er.source();
            let target = er.target();

            // add an edge
            // from "the last node of the source"
            // to "the first node of the target"
            let (_, last_of_source) = m.get(&source).unwrap();
            let (first_of_target, _) = m.get(&target).unwrap();
            graph.add_edge(*last_of_source, *first_of_target, SimpleSeqEdge::new(None));
        }

        SeqGraph(graph)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn genome_graph_linear() {
        let mut g = DiGraph::new();
        g.add_node(GenomeNode::new(b"GATGCTAGTT", 1));
        let gg = GenomeGraph(g);
        let sg = gg.to_seq_graph();
        assert_eq!(sg.node_count(), 10);
        assert_eq!(sg.edge_count(), 9);
    }

    #[test]
    fn genome_graph_branch() {
        let mut g = DiGraph::new();
        let v1 = g.add_node(GenomeNode::new(b"ATCGG", 1));
        let v2 = g.add_node(GenomeNode::new(b"AAC", 1));
        let v3 = g.add_node(GenomeNode::new(b"TTCG", 2));
        g.add_edge(v1, v3, GenomeEdge::new(None));
        g.add_edge(v2, v3, GenomeEdge::new(None));
        let gg = GenomeGraph(g);
        let sg = gg.to_seq_graph();
        assert_eq!(sg.node_count(), 12);
        assert_eq!(sg.edge_count(), 11);
    }

    #[test]
    fn genome_graph_circular() {
        let mut g = DiGraph::new();
        let v1 = g.add_node(GenomeNode::new(b"ATCGG", 1));
        g.add_edge(v1, v1, GenomeEdge::new(None));
        let gg = GenomeGraph(g);
        let sg = gg.to_seq_graph();
        assert_eq!(sg.node_count(), 5);
        assert_eq!(sg.edge_count(), 5);
    }
}
