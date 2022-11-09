//!
//!
//!
//!
use petgraph::{
    graph::{DiGraph, EdgeIndex, Graph, NodeIndex, UnGraph},
    Direction,
};

///
/// remove nodes whose in/out-degree is 1
///
pub fn compact_simple_paths<N: Clone, E: Clone>(graph: &DiGraph<N, E>) -> DiGraph<N, Vec<E>> {
    let mut compacted = graph.map(
        |_node, node_weight| node_weight.clone(),
        |_edge, edge_weight| vec![edge_weight.clone()],
    );
    for node in compacted.node_indices() {
        let in_degree = compacted.edges_directed(node, Direction::Incoming).count();
        let out_degree = compacted.edges_directed(node, Direction::Outgoing).count();
        if in_degree == 1 && out_degree == 1 {
            let in_edge = compacted
                .edges_directed(node, Direction::Incoming)
                .next()
                .unwrap();
            let out_edge = compacted
                .edges_directed(node, Direction::Outgoing)
                .next()
                .unwrap();
            compacted.remove_node(node);
        }
    }
    compacted
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::dot::Dot;

    #[test]
    fn compact_simple_paths_test_1() {
        let g: DiGraph<(), ()> =
            DiGraph::from_edges(&[(0, 1), (1, 2), (2, 4), (0, 3), (3, 4), (4, 0)]);
        let h = compact_simple_paths(&g);
        println!("{:?}", Dot::with_config(&h, &[]));
    }
}
