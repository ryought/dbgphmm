//!
//! `GenomeGraph`
//! Node: linear sequence
//! Edge: its adjacency
//!
use super::seq_graph::{SimpleSeqEdge, SimpleSeqGraph, SimpleSeqNode};
use crate::common::{CopyNum, Reads, Sequence};
use crate::graph::seq_graph::SeqGraph;
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::sample::SampleProfile;
use itertools::Itertools;
use petgraph::dot::Dot;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use petgraph::visit::IntoNodeReferences;
use std::collections::HashMap;

/// Read sampling profile
///
#[derive(Clone, Debug)]
pub struct ReadProfile {
    pub sample_profile: SampleProfile,
    pub phmm_params: PHMMParams,
}

/// GenomeGraph
pub struct GenomeGraph(pub DiGraph<GenomeNode, GenomeEdge>);

/// Node in GenomeGraph is a copy-number-assigned sequence fragment
pub struct GenomeNode {
    /// sequence fragment of this node
    pub seq: Sequence,
    /// copy number of the sequence
    pub copy_num: CopyNum,
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
    pub copy_num: Option<CopyNum>,
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
    ///
    /// Create `GenomeGraph` from the sequences
    ///
    pub fn from_seqs(seqs: &[Sequence]) -> Self {
        let mut g = DiGraph::new();
        for seq in seqs {
            g.add_node(GenomeNode::new(seq, 1));
        }
        GenomeGraph(g)
    }
    /// Get the number of nodes in the GenomeGraph
    pub fn node_count(&self) -> usize {
        self.0.node_count()
    }
    /// Get the number of edges in the GenomeGraph
    pub fn edge_count(&self) -> usize {
        self.0.edge_count()
    }
    ///
    /// Convert `GenomeGraph` into `SimpleSeqGraph`
    /// split a node containing bases into multiple nodes corresponding each bases
    ///
    pub fn to_seq_graph(&self) -> SimpleSeqGraph {
        let mut graph = DiGraph::new();

        // for each node
        let mut m: HashMap<NodeIndex, (NodeIndex, NodeIndex)> = HashMap::new();
        for (node, node_weight) in self.0.node_references() {
            let seq = &node_weight.seq;
            let copy_num = node_weight.copy_num;

            // add a new node for each bases
            // TODO move this to seq_graph.rs?
            let nodes: Vec<NodeIndex> = seq
                .iter()
                .map(|&base| graph.add_node(SimpleSeqNode::new(copy_num, base)))
                .collect();

            // add edges between adjacent two bases
            let edges: Vec<EdgeIndex> = nodes
                .iter()
                .tuple_windows()
                .map(|(&v0, &v1)| graph.add_edge(v0, v1, SimpleSeqEdge::new(Some(copy_num))))
                .collect();

            // store the first/last node
            m.insert(node, (*nodes.first().unwrap(), *nodes.last().unwrap()));
        }

        // for each node
        for er in self.0.edge_references() {
            let source = er.source();
            let target = er.target();
            let weight = er.weight();

            // add an edge
            // from "the last node of the source"
            // to "the first node of the target"
            let (_, last_of_source) = m.get(&source).unwrap();
            let (first_of_target, _) = m.get(&target).unwrap();
            graph.add_edge(
                *last_of_source,
                *first_of_target,
                SimpleSeqEdge::new(weight.copy_num),
            );
        }

        graph
    }
    /// Sample reads from the genome graph.
    pub fn sample_reads(&self, prof: &ReadProfile) -> Reads {
        // convert to phmm
        let phmm = self.to_seq_graph().to_phmm(prof.phmm_params.clone());

        // sample reads using profile
        phmm.sample_reads(&prof.sample_profile)
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

    #[test]
    fn genome_graph_from_seqs() {
        let seqs = vec![b"ATCGATTCGAT".to_vec(), b"CTCTTCTTCTCT".to_vec()];
        let gg = GenomeGraph::from_seqs(&seqs);
        assert_eq!(gg.node_count(), 2);
        assert_eq!(gg.edge_count(), 0);
    }
}
