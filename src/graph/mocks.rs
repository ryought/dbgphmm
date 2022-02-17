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

/// Example small grpah
///
/// ```text
///     a(x2) ----   ------ c(x2)
///               \ /
///                X
///               / \
///     b(x2) ----   ------ d(x2)
/// ```
///
pub fn mock_crossing(has_edge_copy_number: bool) -> GenomeGraph {
    let mut g = DiGraph::new();
    let v1 = g.add_node(GenomeNode::new(&generate(10, 2), 2));
    let v2 = g.add_node(GenomeNode::new(&generate(10, 3), 2));
    let v3 = g.add_node(GenomeNode::new(&generate(10, 4), 2));
    let v4 = g.add_node(GenomeNode::new(&generate(10, 5), 2));
    if has_edge_copy_number {
        g.add_edge(v1, v3, GenomeEdge::new(Some(2)));
        g.add_edge(v2, v3, GenomeEdge::new(Some(0)));
        g.add_edge(v1, v4, GenomeEdge::new(Some(0)));
        g.add_edge(v2, v4, GenomeEdge::new(Some(2)));
    } else {
        g.add_edge(v1, v3, GenomeEdge::new(None));
        g.add_edge(v2, v3, GenomeEdge::new(None));
        g.add_edge(v1, v4, GenomeEdge::new(None));
        g.add_edge(v2, v4, GenomeEdge::new(None));
    }
    GenomeGraph(g)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ni;

    #[test]
    fn test_mock_linear_random() {
        let gg = mock_linear_random(100, 1);
        println!("{}", gg);
        assert_eq!(gg.node_count(), 1);
        assert_eq!(gg.edge_count(), 0);
    }

    #[test]
    fn test_mock_crossing() {
        let g = mock_crossing(true);
        println!("{}", g);
        assert_eq!(g.0[ni(0)].seq, b"TGCTCTGGCG");
        assert_eq!(g.0[ni(1)].seq, b"ATTAGGAGCA");
        assert_eq!(g.0[ni(2)].seq, b"GCTGATAGGG");
        assert_eq!(g.0[ni(3)].seq, b"CGAAGATGAG");
    }
}
