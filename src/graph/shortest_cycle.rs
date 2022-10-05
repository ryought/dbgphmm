//!
//! shortest cycle problem
//!
use petgraph::algo::astar;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;

///
/// calculate the shortest cycle that has an edge `edge_with = (v -> w)` and does not have an edge
/// `edge_without = (w -> v)` (i.e. the reverse edge of `edge_with`) in the cycle.
///
/// current implementation uses astar algorithm
///
pub fn shortest_cycle<N, E>(
    graph: &DiGraph<N, E>,
    edge_with: EdgeIndex,
    edge_without: Option<EdgeIndex>,
) -> Option<Vec<NodeIndex>> {
    let (v, w) = graph.edge_endpoints(edge_with).unwrap();
    if let Some(edge_without) = edge_without {
        let (w2, v2) = graph.edge_endpoints(edge_without).unwrap();
        assert_eq!(v, v2);
        assert_eq!(w, w2);
    }
    // search for shortest path w -> v
    match astar(
        graph,
        w,
        |n| n == v,
        |e| {
            if edge_without.is_some() && e.id() == edge_without.unwrap() {
                // edge_without is not usable
                usize::MAX
            } else {
                // all edges except edge_without have a unit cost
                1
            }
        },
        |_| 0,
    ) {
        Some((cost, path)) => Some(path),
        None => None,
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
    use petgraph::visit::Bfs;

    #[test]
    fn starting() {
        let graph =
            DiGraph::<(), ()>::from_edges(&[(0, 1), (1, 0), (1, 2), (2, 1), (2, 0), (0, 2)]);
        let cycle = shortest_cycle(&graph, EdgeIndex::new(0), Some(EdgeIndex::new(1)));
        println!("{:?}", cycle);
        assert_eq!(
            cycle,
            Some(vec![
                NodeIndex::new(1),
                NodeIndex::new(2),
                NodeIndex::new(0)
            ])
        )
    }
}
