use super::genome_graph::{GenomeEdge, GenomeGraph, GenomeNode};
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

pub fn mock() -> GenomeGraph {
    let mut g = DiGraph::new();
    let v1 = g.add_node(GenomeNode::new(b"ATTCGATCGTCG", 1));
    let v2 = g.add_node(GenomeNode::new(b"ATCGATG", 1));
    let v3 = g.add_node(GenomeNode::new(b"TTCGAT", 2));
    g.add_edge(v1, v3, GenomeEdge::new(None));
    g.add_edge(v2, v3, GenomeEdge::new(None));
    GenomeGraph(g)
}
