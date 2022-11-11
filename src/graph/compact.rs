//!
//! Graph compacting functions
//! Functions in this file handles generic petgraph::DiGraph, not Dbg struct.
//!
//! * remove_deadends
//! * collapse_simple_paths
//!
use petgraph::{
    graph::{DiGraph, EdgeIndex, Graph, NodeIndex, UnGraph},
    Direction,
};

///
/// get node degree (in degree if dir=Incoming, out degree if dir=Outgoing)
///
fn node_degree<N, E>(graph: &DiGraph<N, E>, node: NodeIndex, dir: Direction) -> usize {
    graph.edges_directed(node, dir).count()
}

///
/// find a deadend node in graph (in order of index)
///
fn find_deadend<N, E>(graph: &DiGraph<N, E>) -> Option<NodeIndex> {
    graph.node_indices().find(|&node| {
        node_degree(graph, node, Direction::Incoming) == 0
            || node_degree(graph, node, Direction::Outgoing) == 0
    })
}

///
/// Remove all deadend nodes
/// * in_degree is zero
/// * out_degree is zero
///
/// ## Notes
/// * Removing can cause node index change
/// * There is room to improve, current algorithm takes O(N^2).
///
pub fn remove_deadends<N, E>(mut graph: DiGraph<N, E>) -> DiGraph<N, E> {
    while let Some(deadend) = find_deadend(&graph) {
        graph.remove_node(deadend);
    }
    graph
}

///
/// WIP
///
/// remove nodes whose in/out-degree is 1
///
pub fn compact_simple_paths<N: Clone, E: Clone>(
    graph: &DiGraph<N, E>,
) -> DiGraph<N, Vec<(EdgeIndex, E)>> {
    let mut compacted = graph.map(
        |_node, node_weight| node_weight.clone(),
        |edge, edge_weight| vec![(edge, edge_weight.clone())],
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
    use crate::common::{ei, ni};
    use petgraph::dot::Dot;

    #[test]
    fn compact_simple_paths_test_1() {
        let g: DiGraph<(), ()> =
            DiGraph::from_edges(&[(0, 1), (1, 2), (2, 4), (0, 3), (3, 4), (4, 0)]);
        println!("{:?}", Dot::with_config(&g, &[]));
        let h = compact_simple_paths(&g);
        println!("{:?}", Dot::with_config(&h, &[]));
    }

    #[test]
    fn remove_deadends_test() {
        // graph with deadends
        {
            let g: DiGraph<(), ()> = DiGraph::from_edges(&[
                // main circle <1,2,3,4>
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 1),
                //deadend edges 0,5,6
                (2, 0),
                (4, 5),
                (6, 5),
            ]);
            println!("{:?}", Dot::with_config(&g, &[]));
            assert_eq!(node_degree(&g, ni(4), Direction::Incoming), 1);
            assert_eq!(node_degree(&g, ni(4), Direction::Outgoing), 2);
            assert_eq!(node_degree(&g, ni(0), Direction::Incoming), 1);
            assert_eq!(node_degree(&g, ni(0), Direction::Outgoing), 0);
            assert_eq!(find_deadend(&g), Some(ni(0)));

            let h = remove_deadends(g.clone());
            println!("{:?}", Dot::with_config(&h, &[]));
            assert_eq!(h.node_count(), 4);
            assert_eq!(h.edge_count(), 4);
        }

        // graph with deadends
        {
            let g: DiGraph<(), ()> = DiGraph::from_edges(&[(0, 1), (1, 2), (2, 3), (3, 0)]);
            println!("{:?}", Dot::with_config(&g, &[]));
            assert_eq!(find_deadend(&g), None);

            let h = remove_deadends(g.clone());
            println!("{:?}", Dot::with_config(&h, &[]));
            assert_eq!(h.node_count(), 4);
            assert_eq!(h.edge_count(), 4);
        }
    }
}
