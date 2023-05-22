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
/// Get 1st to k-th shortest cycle starting from `edge`
///
/// [k shortest path routing - wikipedia](https://en.wikipedia.org/wiki/K_shortest_path_routing)
///
///
/// # Algorithm
///
/// `heap`: Candidate set of shortest paths starting from `edge` sorted by increasing order of
/// scores.
///
/// `count[v]`: How many times the accepted shortest path visited the node `v`
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
    let (source, _) = graph.edge_endpoints(edge).unwrap();
    let mut cycles = vec![];

    while let Some(PathWithScore(score, path)) = heap.pop() {
        //
        // This `path` is shortest path from source to node u with `score`.
        //
        let last_edge = *path.last().unwrap();
        let (_, u) = graph.edge_endpoints(last_edge).unwrap();
        //
        // mark as visited
        //
        count[u.index()] += 1;
        if u == source && count[u.index()] <= k {
            cycles.push(path.clone());
        }
        //
        // if u is visited more than k times, extending is redundant
        // (only <=k shortest path is used in k-th shortest path from target to source)
        //
        if count[u.index()] <= k {
            for e in graph.edges_directed(u, Direction::Outgoing) {
                if is_joinable(&path, e.id()) {
                    let mut path = path.clone();
                    path.push(e.id());
                    heap.push(PathWithScore(score + edge_cost(e.id()), path))
                }
            }
        }
    }
    // println!("{} cycles found", cycles.len());

    cycles
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ei;
    use rustflow::min_flow::residue::{
        is_meaningful_move_on_residue_graph, ResidueDirection, ResidueEdge, ResidueGraph,
    };

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

        let r = k_shortest_cycle(&graph, ei(0), 1, |_| 1, |_, _| true);
        println!("{:?}", r);
        println!("{:?}", r.len());
        assert_eq!(r, vec![vec![ei(0), ei(1), ei(2), ei(3)]]);

        let r = k_shortest_cycle(&graph, ei(4), 3, |_| 1, |_, _| true);
        println!("{:?}", r);
        println!("{:?}", r.len());
        assert_eq!(
            r,
            vec![
                vec![ei(4), ei(6), ei(7), ei(8)],
                vec![ei(4), ei(6), ei(7), ei(8), ei(4), ei(6), ei(7), ei(8)],
                vec![
                    ei(4),
                    ei(6),
                    ei(7),
                    ei(8),
                    ei(4),
                    ei(6),
                    ei(7),
                    ei(8),
                    ei(4),
                    ei(6),
                    ei(7),
                    ei(8)
                ],
            ]
        );
    }
    #[test]
    fn k_shortest_02() {
        let up = ResidueDirection::Up;
        let down = ResidueDirection::Down;
        let rg: ResidueGraph<usize> = Graph::from_edges(&[
            (2, 5, ResidueEdge::new(1, 0.0, ei(0), up)),
            (3, 5, ResidueEdge::new(1, 0.0, ei(1), up)),
            (5, 3, ResidueEdge::new(1, 0.0, ei(1), down)),
            (6, 3, ResidueEdge::new(1, 0.0, ei(2), up)),
            (3, 6, ResidueEdge::new(1, 0.0, ei(2), down)),
            (2, 4, ResidueEdge::new(1, 0.0, ei(3), up)),
            (4, 2, ResidueEdge::new(1, 0.0, ei(3), down)),
            (3, 4, ResidueEdge::new(1, 0.0, ei(4), up)),
            (0, 6, ResidueEdge::new(1, 0.0, ei(5), up)),
            (6, 0, ResidueEdge::new(1, 0.0, ei(5), down)),
            (4, 1, ResidueEdge::new(1, 0.0, ei(6), up)),
            (1, 4, ResidueEdge::new(1, 0.0, ei(6), down)),
            (1, 0, ResidueEdge::new(1, 0.0, ei(7), up)),
            (0, 1, ResidueEdge::new(1, 0.0, ei(7), down)),
            (6, 2, ResidueEdge::new(1, 0.0, ei(8), up)),
            (2, 6, ResidueEdge::new(1, 0.0, ei(8), down)),
            (5, 1, ResidueEdge::new(1, 0.0, ei(9), up)),
            (1, 5, ResidueEdge::new(1, 0.0, ei(9), down)),
        ]);
        let weights = vec![1, 1, 5170, 1, 1, 403, 331, 49, 5158, 331];

        // Mimics compacted DBG for finding rescue cycles
        // Residue graph
        //     2 -> 5 [ label = "ResidueEdge { count: 1, weight: 0.0, target: EdgeIndex(0), direction: Up }" ]
        //     3 -> 5 [ label = "ResidueEdge { count: 1, weight: 0.0, target: EdgeIndex(1), direction: Up }" ]
        //     5 -> 3 [ label = "ResidueEdge { count: 1, weight: -0.0, target: EdgeIndex(1), direction: Down }" ]
        //     6 -> 3 [ label = "ResidueEdge { count: 1, weight: 0.0, target: EdgeIndex(2), direction: Up }" ]
        //     3 -> 6 [ label = "ResidueEdge { count: 1, weight: -0.0, target: EdgeIndex(2), direction: Down }" ]
        //     2 -> 4 [ label = "ResidueEdge { count: 1, weight: 0.0, target: EdgeIndex(3), direction: Up }" ]
        //     4 -> 2 [ label = "ResidueEdge { count: 1, weight: -0.0, target: EdgeIndex(3), direction: Down }" ]
        //     3 -> 4 [ label = "ResidueEdge { count: 1, weight: 0.0, target: EdgeIndex(4), direction: Up }" ]
        //     0 -> 6 [ label = "ResidueEdge { count: 1, weight: 0.0, target: EdgeIndex(5), direction: Up }" ]
        //     6 -> 0 [ label = "ResidueEdge { count: 1, weight: -0.0, target: EdgeIndex(5), direction: Down }" ]
        //     4 -> 1 [ label = "ResidueEdge { count: 1, weight: 0.0, target: EdgeIndex(6), direction: Up }" ]
        //     1 -> 4 [ label = "ResidueEdge { count: 1, weight: -0.0, target: EdgeIndex(6), direction: Down }" ]
        //     1 -> 0 [ label = "ResidueEdge { count: 1, weight: 0.0, target: EdgeIndex(7), direction: Up }" ]
        //     0 -> 1 [ label = "ResidueEdge { count: 1, weight: -0.0, target: EdgeIndex(7), direction: Down }" ]
        //     6 -> 2 [ label = "ResidueEdge { count: 1, weight: 0.0, target: EdgeIndex(8), direction: Up }" ]
        //     2 -> 6 [ label = "ResidueEdge { count: 1, weight: -0.0, target: EdgeIndex(8), direction: Down }" ]
        //     5 -> 1 [ label = "ResidueEdge { count: 1, weight: 0.0, target: EdgeIndex(9), direction: Up }" ]
        //     1 -> 5 [ label = "ResidueEdge { count: 1, weight: -0.0, target: EdgeIndex(9), direction: Down }" ]
        //
        // with weights
        //     e0 1
        //     e1 1
        //     e2 5170
        //     e3 1
        //     e4 1
        //     e5 403
        //     e6 331
        //     e7 49
        //     e8 5158
        //     e9 331
        //
        // e7,e10
        let k = 10;
        let cycles = k_shortest_cycle(
            &rg,
            EdgeIndex::new(7),
            k,
            |e| weights[rg[e].target.index()],
            |path, edge| {
                let last_edge = path.last().copied().unwrap();
                is_meaningful_move_on_residue_graph(&rg, last_edge, edge)
            },
        );

        for (i, cycle) in cycles.iter().enumerate() {
            println!("{} {:?}", i, cycle);
        }
    }
}
