use super::genome_graph::{GenomeEdge, GenomeGraph, GenomeNode};
use crate::random_seq::generate;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

///
/// Sample linear genome (10bp) "ATTCGATCGT"
///
pub fn mock_linear() -> GenomeGraph {
    let mut g = DiGraph::new();
    g.add_node(GenomeNode::new(b"ATTCGATCGT", 1));
    GenomeGraph(g)
}

///
/// Random linear genome with length and from seed
///
pub fn mock_linear_random(length: usize, seed: u64) -> GenomeGraph {
    let seq = generate(length, seed);
    let mut g = DiGraph::new();
    g.add_node(GenomeNode::new(&seq, 1));
    GenomeGraph(g)
}

pub fn mock() -> GenomeGraph {
    let mut g = DiGraph::new();
    let v1 = g.add_node(GenomeNode::new(b"ATTCGATCGTCG", 1));
    let v2 = g.add_node(GenomeNode::new(b"ATCGATG", 1));
    let v3 = g.add_node(GenomeNode::new(b"TTCGAT", 2));
    g.add_edge(v1, v3, GenomeEdge::new(None));
    g.add_edge(v2, v3, GenomeEdge::new(None));
    GenomeGraph(g)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mock_linear_random() {
        let gg = mock_linear_random(100, 1);
        println!("{}", gg);
        assert_eq!(gg.node_count(), 1);
        assert_eq!(gg.edge_count(), 0);
    }
}
