//!
//! k shortest path algorithm
//!
//! petgraph::algo::k_shortest_path::k_shortest_path
//!
use petgraph::{
    graph::{DiGraph, EdgeIndex, Graph, NodeIndex, UnGraph},
    visit::EdgeRef,
    Direction,
};
use std::collections::BinaryHeap;

//
// Internal Path with score struct
//

#[derive(Clone, Debug)]
struct PathWithScore(usize, Vec<EdgeIndex>);

// reverse order
impl Ord for PathWithScore {
    fn cmp(&self, other: &PathWithScore) -> std::cmp::Ordering {
        let a = &self.0;
        let b = &other.0;
        b.cmp(a)
    }
}

impl PartialEq for PathWithScore {
    fn eq(&self, other: &PathWithScore) -> bool {
        self.cmp(other) == std::cmp::Ordering::Equal
    }
}

impl Eq for PathWithScore {}

impl PartialOrd for PathWithScore {
    fn partial_cmp(&self, other: &PathWithScore) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

//
//
//

///
/// Get k-th shortest cycle
///
/// [k shortest path routing - wikipedia](https://en.wikipedia.org/wiki/K_shortest_path_routing)
///
pub fn k_shortest_cycle<N, E, F, G>(
    graph: &DiGraph<N, E>,
    edge: EdgeIndex,
    k: usize,
    edge_cost: F,
    is_joinable: G,
) -> Vec<Vec<EdgeIndex>>
where
    F: Fn(EdgeIndex) -> usize,
    G: Fn(&[EdgeIndex], EdgeIndex) -> bool,
{
    let mut count: Vec<usize> = vec![0; graph.node_count()];
    let mut heap = BinaryHeap::new();
    heap.push(PathWithScore(edge_cost(edge), vec![edge]));
    let (s, _) = graph.edge_endpoints(edge).unwrap();
    let mut cycles = vec![];

    while let Some(PathWithScore(score, path)) = heap.pop() {
        let (_, u) = graph.edge_endpoints(*path.last().unwrap()).unwrap();
        count[u.index()] += 1;
        if u == s {
            cycles.push(path.clone());
        }
        if count[u.index()] < k {
            for e in graph.edges_directed(u, Direction::Outgoing) {
                if is_joinable(&path, e.id()) {
                    let mut path = path.clone();
                    path.push(e.id());
                    heap.push(PathWithScore(score + edge_cost(e.id()), path))
                }
            }
        }
    }

    cycles
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ei;

    #[test]
    fn path_with_score() {
        let mut h = BinaryHeap::new();
        let a = PathWithScore(10, vec![ei(5), ei(1)]);
        let b = PathWithScore(3, vec![ei(4), ei(2)]);
        let c = PathWithScore(15, vec![ei(3)]);
        h.push(a);
        h.push(b);
        h.push(c);
        while let Some(p) = h.pop() {
            println!("p={:?}", p);
        }
    }
    #[test]
    fn k_shortest_01() {
        let mut graph: DiGraph<(), ()> = Graph::new();
        let a = graph.add_node(()); // node with no weight
        let b = graph.add_node(());
        let c = graph.add_node(());
        let d = graph.add_node(());
        let e = graph.add_node(());
        let f = graph.add_node(());
        let g = graph.add_node(());
        let h = graph.add_node(());
        // z will be in another connected component
        let z = graph.add_node(());

        graph.extend_with_edges(&[
            (a, b), //e0
            (b, c),
            (c, d),
            (d, a),
            (e, f),
            (b, e), //e5
            (f, g),
            (g, h),
            (h, e),
        ]);
        // a ----> b ----> e ----> f
        // ^       |       ^       |
        // |       v       |       v
        // d <---- c       h <---- g

        let r = k_shortest_cycle(&graph, ei(0), 4, |_| 1, |_, _| true);
        println!("{:?}", r);
    }
}
