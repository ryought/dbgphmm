//!
//! k shortest path algorithm
//!
//! petgraph::algo::k_shortest_path::k_shortest_path
//!
use petgraph::{
    graph::{DiGraph, EdgeIndex, Graph, NodeIndex, UnGraph},
    stable_graph::StableDiGraph,
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

/// Find k-shortest paths from source to target without edge-repetition (i.e. simple path)
///
/// by using `petgraph::algo::astar::astar`
///
///
/// # TODO
///
/// * add `is_joinable` condition argument?
///
pub fn k_shortest_simple_path<N: Clone, E: Clone, F>(
    graph: &DiGraph<N, E>,
    source: NodeIndex,
    target: NodeIndex,
    k: usize,
    edge_cost: F,
    // is_joinable: G,
) -> Vec<Vec<EdgeIndex>>
where
    F: Fn(EdgeIndex) -> usize,
    // G: Fn(&[EdgeIndex], EdgeIndex) -> bool,
{
    // List A
    let mut paths: Vec<Vec<EdgeIndex>> = Vec::new();
    // List B
    let mut candidates: BinaryHeap<PathWithScore> = BinaryHeap::new();
    // let mut candidates: Vec<PathWithScore> = Vec::new();

    // Iteration 0
    let graph = StableDiGraph::from(graph.clone());
    let (_, a0) = shortest_path(&graph, source, target, &edge_cost)
        .expect("no path between source and target");
    paths.push(a0);

    // Iteration k
    for iter in 1..k {
        // (i) create candidates of a_k from a_0 ... a_k-1
        let mut graph = graph.clone();
        let a = &paths[paths.len() - 1];

        // println!("k={} a[k-1]={:?}", iter, a);
        // for each edge a[i] in path a
        // create a candidate path a[0]..a[i-1] + different edge from a[i]
        for i in 0..a.len() {
            // a. create graph with new weight
            let (v, _) = graph.edge_endpoints(a[i]).unwrap();

            // (1) delete i to i+1
            for aj in paths.iter() {
                if aj.len() >= i && aj[..i] == a[..i] {
                    // println!("removing edge {:?}", aj[i]);
                    graph.remove_edge(aj[i]);
                }
            }
            // b. find shortest path from v to terminal and store it in List B
            match shortest_path(&graph, v, target, &edge_cost) {
                Some((_, path)) => {
                    let mut ak = Vec::new();
                    ak.extend_from_slice(&a[..i]); // path from s to v
                    ak.extend_from_slice(&path); // path from v to t
                    let score = total_cost(&graph, &ak, &edge_cost);
                    // println!("candidate {:?} {} ", ak, score);
                    if candidates.iter().all(|PathWithScore(_, path)| path != &ak) {
                        // println!("new");
                        candidates.push(PathWithScore(score, ak));
                    } else {
                        // println!("old");
                    }
                }
                None => {
                    // println!("no candidate");
                }
            }

            // (2) delete node 1 to i-1
            // println!("removing node {:?}", v);
            graph.remove_node(v);
        }

        // (ii)
        // a. pick minimum cost path from List B and add to List A
        match candidates.pop() {
            Some(PathWithScore(_, ak)) => {
                // println!("path {:?}", ak);
                paths.push(ak);
            }
            None => {
                // println!("no more candidates");
                break;
            }
        }
    }

    paths
}

///
/// Compute shortest path from source to target using a-star algorithm in petgraph
///
pub fn shortest_path<N, E, F>(
    graph: &StableDiGraph<N, E>,
    source: NodeIndex,
    target: NodeIndex,
    edge_cost: F,
) -> Option<(usize, Vec<EdgeIndex>)>
where
    F: Fn(EdgeIndex) -> usize,
{
    petgraph::algo::astar::astar(graph, source, |v| v == target, |e| edge_cost(e.id()), |_| 0)
        .map(|(cost, nodes)| (cost, nodes_to_edges(graph, &nodes, edge_cost)))
}

///
/// convert shortest path v -> w as node sequence into edge sequence by picking minimum cost
/// edge between two adjacent nodes.
///
pub fn nodes_to_edges<N, E, F>(
    graph: &StableDiGraph<N, E>,
    nodes: &[NodeIndex],
    edge_cost: F,
) -> Vec<EdgeIndex>
where
    F: Fn(EdgeIndex) -> usize,
{
    let mut edges = Vec::new();
    let n = nodes.len();

    // convert (nodes[i], nodes[i+1]) into an edge
    for i in 0..(n - 1) {
        let v = nodes[i];
        let w = nodes[i + 1];

        // pick a minimum cost edge between v and w
        let edge = graph
            .edges_connecting(v, w)
            .min_by_key(|e| edge_cost(e.id()))
            .unwrap();

        edges.push(edge.id());
    }

    edges
}

///
///
///
pub fn total_cost<N, E, F>(graph: &StableDiGraph<N, E>, edges: &[EdgeIndex], edge_cost: F) -> usize
where
    F: Fn(EdgeIndex) -> usize,
{
    edges.iter().map(|&e| edge_cost(e)).sum()
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{ei, ni};
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

        // let cycles =
        //     k_shortest_simple_path(&rg, EdgeIndex::new(7), k, |e| weights[rg[e].target.index()]);
    }
    #[test]
    fn k_shortest_simple_01() {
        let graph: DiGraph<(), ()> = Graph::from_edges(&[
            (0, 1),
            (1, 2),
            (2, 3),
            // alt of 1->2
            (1, 4),
            (4, 5),
            (5, 2),
            // alt of 2->3
            (2, 6),
            (6, 3),
        ]);
        let paths = k_shortest_simple_path(&graph, ni(0), ni(3), 10, |_| 1);
        println!("paths={:?}", paths);
        assert_eq!(
            paths,
            vec![
                vec![ei(0), ei(1), ei(2)],
                vec![ei(0), ei(1), ei(6), ei(7)],
                vec![ei(0), ei(3), ei(4), ei(5), ei(2)],
                vec![ei(0), ei(3), ei(4), ei(5), ei(6), ei(7)],
            ]
        );
    }
}
