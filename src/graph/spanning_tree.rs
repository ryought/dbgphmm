//!
//! Generate spanning tree by breadth first search from a specified node
//!
use petgraph::graph::{DiGraph, EdgeIndex, Graph, NodeIndex, UnGraph};
use petgraph::visit::{Bfs, VisitMap, Visitable};
use std::collections::VecDeque;

struct EdgeInfo {
    origin: EdgeIndex,
}

struct SpanningTree {
    tree: UnGraph<(), EdgeInfo>,
    edges: Vec<EdgeIndex>,
}

///
/// create a spanning tree of Graph (direction of edges is ignored) by breadth first search from a
/// specified node
///
fn spanning_tree<N, E>(graph: &UnGraph<N, E>, start: NodeIndex) -> SpanningTree {
    // do bfs
    let mut visited = graph.visit_map();
    visited.visit(start);
    let mut queue = VecDeque::new();
    queue.push_back(start);

    while let Some(node) = queue.pop_front() {
        // for edge in graph.edges(node) {
        for succ in graph.neighbors(node) {
            // let graph.edge_endpoints()
            if visited.visit(succ) {
                queue.push_back(succ);
            }
        }
    }

    unimplemented!();
}

///
/// convert directed graph into undirected graph
///
fn to_undirected_graph<N, E>(digraph: DiGraph<N, E>) -> UnGraph<N, E> {
    digraph.into_edge_type()
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::dot::Dot;

    #[test]
    fn spanning_tree_mock_undirected() {
        let mut g: DiGraph<(), ()> = DiGraph::new();
        g.extend_with_edges(&[(0, 1), (2, 0), (2, 1)]);

        let mut g2 = to_undirected_graph(g);
        println!("{:?}", Dot::with_config(&g2, &[]));

        // spanning_tree(&g2, NodeIndex::new(0));
    }
}
