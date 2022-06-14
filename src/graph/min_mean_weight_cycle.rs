//!
//! Find a minimum mean-weight cycle in petgraph::DiGraph
//!
//! Karp's minimum mean-weight cycle algorithm
//!
//! # References
//!
//! * https://walkccc.me/CLRS/Chap24/Problems/24-5/
//!
use super::FloatWeight;
use petgraph::prelude::*;
use petgraph::visit::{VisitMap, Visitable};

//
// basic shortest paths with node indexing
//

///
/// F[k][v] = (min weight path (with length k) from source to a node v)
/// for k=0..n and all nodes v (reachable from source node) n=|V|
///
/// and backtracking information
///
#[derive(Clone, Debug)]
pub struct ShortestPaths {
    ///
    /// Distances
    ///
    /// `dists[k: length of path][v: node]`
    ///
    dists: Vec<Vec<f64>>,
    ///
    /// Predecessors for backtracking
    /// `preds[k: length of path][v: node] = (w: node, e: edge)`
    /// means that "min weight path `F[k][v]` ends with the `e=w->v` edge"
    ///
    preds: Vec<Vec<Option<(NodeIndex, EdgeIndex)>>>,
}

impl ShortestPaths {}

///
/// # Inputs
///
/// Directed graph with edge attribute E has f64 weight.
///
pub fn shortest_paths<N, E: FloatWeight>(
    graph: &DiGraph<N, E>,
    source: NodeIndex,
) -> ShortestPaths {
    let n = graph.node_count();
    let ix = |node: NodeIndex| node.index();

    // (1) Initialize
    let mut dists = vec![vec![f64::INFINITY; n]; n + 1];
    let mut preds = vec![vec![None; n]; n + 1];
    dists[0][ix(source)] = 0.0;

    // (2) Update
    for k in 1..=n {
        for u in graph.node_indices() {
            for edge in graph.edges(u) {
                // edge e: u->v with weight w
                let e = edge.id();
                let v = edge.target();
                let w = edge.weight().float_weight();
                // if s->u->v is shorter than s->v, update the route.
                if dists[k - 1][ix(u)] + w < dists[k][ix(v)] {
                    dists[k][ix(v)] = dists[k - 1][ix(u)] + w;
                    preds[k][ix(v)] = Some((u, e));
                }
            }
        }
    }

    ShortestPaths { dists, preds }
}

//
// custom shortest paths with edge indexing for using edge-adjacency condition
//

///
/// F[k][e] = (min weight path s.t.
///   * length is `k`
///   * ends with edge `e`
///   * from source to a target node of `e`)
///
/// for `k=1..n` and all edges `e` (reachable from source node)
/// (Here `n=|V|` because path length >= `n` means node repetition, so it should have a cycle)
///
/// and backtracking information
///
#[derive(Clone, Debug)]
pub struct ShortestPathsByEdge {
    ///
    /// Distances
    ///
    /// `dists[k: length of path][e: edge]`
    ///
    dists: Vec<Vec<f64>>,
    ///
    /// Predecessors for backtracking
    /// `preds[k: length of path][e: edge] = e': edge`
    /// means that "min weight path `F[k][e]` ends with edges `(...,e',e)`"
    ///
    preds: Vec<Vec<Option<EdgeIndex>>>,
}

impl ShortestPathsByEdge {
    ///
    /// Convert edge-indexed ShortestPathsByEdge into
    /// node-indexed ShortestPaths.
    ///
    pub fn into_shortest_paths<N, E>(self, graph: &DiGraph<N, E>) -> ShortestPaths {
        let n = graph.node_count();
        let mut dists = vec![vec![f64::INFINITY; n]; n + 1];
        let mut preds = vec![vec![None; n]; n + 1];

        // TODO add assertion of ShortestPathsByEdge
        assert_eq!(self.dists.len(), n);
        assert_eq!(self.preds.len(), n);
        for k in 0..n {
            // TODO
        }
        ShortestPaths { dists, preds }
    }
}

///
/// # Inputs
///
/// Directed graph with edge attribute E has f64 weight.
///
pub fn shortest_paths_by_edge<N, E, F>(
    graph: &DiGraph<N, E>,
    source: NodeIndex,
    edge_moveable: F,
) -> ShortestPathsByEdge
where
    E: FloatWeight,
    F: Fn(EdgeIndex, EdgeIndex) -> bool,
{
    let n_nodes = graph.node_count();
    let n_edges = graph.edge_count();
    let ix = |edge: EdgeIndex| edge.index();
    let eps = f64::EPSILON;

    // (1) Initialize
    //     e
    // s ----> *
    let mut dists = vec![vec![f64::INFINITY; n_edges]; n_nodes];
    let mut preds = vec![vec![None; n_edges]; n_nodes];
    for edge in graph.edges_directed(source, Direction::Outgoing) {
        dists[0][ix(edge.id())] = edge.weight().float_weight();
    }

    // (2) Update
    for k in 1..n_nodes {
        // for each edge
        // * weight w
        // * from a to b
        for edge in graph.edge_references() {
            let v = edge.source();
            let w = edge.weight().float_weight();
            let e = edge.id();

            for edge_pred in graph.edges_directed(v, Direction::Incoming) {
                let ep = edge_pred.id();

                if edge_moveable(ep, e) {
                    //  ep       e
                    // ----> v ----> *
                    if dists[k - 1][ix(ep)] + w + eps < dists[k][ix(e)] {
                        dists[k][ix(e)] = dists[k - 1][ix(ep)] + w;
                        preds[k][ix(e)] = Some((ep));
                    }
                }
            }
        }
    }

    ShortestPathsByEdge { dists, preds }
}

///
/// Find a minimizer pair `(k, v)`
/// that satisfies `min_v max_k (F[n][v] - F[k][v]) / (n-k)`.
///
pub fn find_minimizer_pair<N, E: FloatWeight>(
    graph: &DiGraph<N, E>,
    source: NodeIndex,
    paths: &ShortestPaths,
) -> Option<(usize, NodeIndex, f64)> {
    let n = graph.node_count();
    (0..n)
        .filter_map(|v| {
            (0..n)
                .filter_map(|k| {
                    let fnv = paths.dists[n][v];
                    let fkv = paths.dists[k][v];
                    if fnv != f64::INFINITY && fkv != f64::INFINITY {
                        Some((k, (fnv - fkv) / (n - k) as f64))
                    } else {
                        None
                    }
                })
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                .map(|(k, score)| (k, NodeIndex::new(v), score))
        })
        .min_by(|(_, _, a), (_, _, b)| a.partial_cmp(b).unwrap())
}

///
/// Traceback a path from source to target
/// by using `preds` in paths
///
pub fn traceback_preds<N, E: FloatWeight>(
    graph: &DiGraph<N, E>,
    source: NodeIndex,
    target: NodeIndex,
    paths: &ShortestPaths,
) -> (Vec<NodeIndex>, Vec<EdgeIndex>) {
    let n = graph.node_count();
    let ix = |node: NodeIndex| node.index();

    let mut node = target;
    let mut nodes = vec![target];
    let mut edges = vec![];

    for k in (1..=n).rev() {
        //            edge_pred
        // node_pred -----------> node
        let (node_pred, edge_pred) = match paths.preds[k][ix(node)] {
            Some((node_pred, edge_pred)) => (node_pred, edge_pred),
            None => panic!("no parent"),
        };
        nodes.push(node_pred);
        edges.push(edge_pred);
        node = node_pred;
    }

    nodes.reverse();
    edges.reverse();
    (nodes, edges)
}

///
/// Find a cycle in a path
///
pub fn find_cycle<N, E>(graph: &DiGraph<N, E>, path: &[NodeIndex]) -> Option<Vec<NodeIndex>> {
    let mut v = graph.visit_map();
    let mut cycle = Vec::new();

    for &node in path.iter().rev() {
        if v.is_visited(&node) {
            // this node is visited twice.  cycle is detected.
            let prev_visit = cycle.iter().position(|&v| v == node).unwrap();
            cycle = cycle[prev_visit..].to_vec();
            cycle.reverse();
            return Some(cycle);
        } else {
            // fisrt visit of this node
            v.visit(node);
            cycle.push(node);
        }
    }

    None
}

///
/// Find a minimum mean-weight cycle in a graph
/// Returns None if there is no cycle.
///
pub fn find_minimum_mean_weight_cycle<N, E: FloatWeight>(
    graph: &DiGraph<N, E>,
    source: NodeIndex,
) -> Option<(Vec<NodeIndex>, f64)> {
    let sp = shortest_paths(graph, source);
    match find_minimizer_pair(graph, source, &sp) {
        Some((_, v, mean_weight)) => {
            let (path, _) = traceback_preds(graph, source, v, &sp);
            let cycle = find_cycle(graph, &path)
                .expect("minimizer pair was found, but no cycle was found when tracebacking");
            Some((cycle, mean_weight))
        }
        None => None,
    }
}

///
/// Find a negative cycle by using `find_minimum_mean_weight_cycle`.
///
pub fn find_negative_cycle<N, E: FloatWeight>(
    graph: &DiGraph<N, E>,
    source: NodeIndex,
) -> Option<Vec<NodeIndex>> {
    match find_minimum_mean_weight_cycle(graph, source) {
        Some((cycle, mean_weight)) => {
            if mean_weight < 0.0 {
                Some(cycle)
            } else {
                None
            }
        }
        None => None,
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{ei, ni};

    #[test]
    fn shortest_paths_00() {
        let mut g: DiGraph<(), f64> = DiGraph::new();
        g.extend_with_edges(&[
            (0, 1, 1.0),
            (0, 2, 1.0),
            (1, 3, 2.0),
            (2, 3, 1.0),
            (2, 4, 2.0),
        ]);
        let sp = shortest_paths(&g, ni(0));
        println!("{:?}", sp);
        assert_eq!(sp.dists[0][0], 0.0);
        assert_eq!(sp.preds[0][0], None);
        // 0->1
        assert_eq!(sp.dists[1][1], 1.0);
        assert_eq!(sp.preds[1][1], Some((ni(0), ei(0))));
        // 0->2
        assert_eq!(sp.dists[1][2], 1.0);
        assert_eq!(sp.preds[1][2], Some((ni(0), ei(1))));
        // 0->2->3
        assert_eq!(sp.dists[2][3], 1.0 + 1.0);
        assert_eq!(sp.preds[2][3], Some((ni(2), ei(3))));
        // 0->2->4
        assert_eq!(sp.dists[2][4], 1.0 + 2.0);
        assert_eq!(sp.preds[2][4], Some((ni(2), ei(4))));

        let cycle = find_minimum_mean_weight_cycle(&g, ni(0));
        println!("cycle={:?}", cycle);
        assert_eq!(cycle, None);
    }

    /// more simple graph with two cycles
    ///         3 ----> 4
    ///         ^       |
    ///         |       V
    /// 0 ----> 1 <---- 5
    /// ^       |
    /// +-- 2 <-+
    #[test]
    fn shortest_paths_01() {
        let mut g: DiGraph<(), f64> = DiGraph::new();
        g.extend_with_edges(&[
            (0, 1, 1.0),
            (1, 2, 3.0),
            (2, 0, 1.0),
            (1, 3, 1.0),
            (3, 4, 2.0),
            (4, 5, 1.0),
            (5, 1, 1.0),
        ]);
        let sp = shortest_paths(&g, ni(0));
        println!("{:?}", sp);
        let (k, v, s) = find_minimizer_pair(&g, ni(0), &sp).unwrap();
        println!("k={} v={} s={}", k, v.index(), s);
        let (path, _) = traceback_preds(&g, ni(0), v, &sp);
        println!("path={:?}", path);
        let cycle = find_cycle(&g, &path);
        println!("cycle={:?}", cycle);

        let cycle = find_minimum_mean_weight_cycle(&g, ni(0));
        println!("cycle={:?}", cycle);
        assert_eq!(cycle, Some((vec![ni(3), ni(4), ni(5), ni(1)], 1.25)));
    }

    /// example graph with parallel edges
    ///     +----->
    /// +-> 0 ----> 1 -----> 3
    /// |           |        |
    /// +----- 2 <--+        |
    ///        ^             |
    ///        +-------------+
    #[test]
    fn shortest_paths_02() {
        let mut g: DiGraph<(), f64> = DiGraph::new();
        g.extend_with_edges(&[
            (0, 1, 3.0),
            (0, 1, 1.0),
            (1, 2, 1.0),
            (2, 0, 1.0),
            (1, 3, 1.0),
            (3, 2, 4.0),
        ]);

        let cycle = find_minimum_mean_weight_cycle(&g, ni(0));
        println!("cycle={:?}", cycle);
        assert_eq!(cycle, Some((vec![ni(2), ni(0), ni(1)], 1.0)));
    }

    ///
    /// counterexample in 'A note on finding minimum mean cycle' (Chaturvedi, 2017)
    ///
    #[test]
    fn shortest_paths_03() {
        let mut g: DiGraph<(), f64> = DiGraph::new();
        g.extend_with_edges(&[
            (0, 1, 1.0),
            (1, 2, 3.0),
            (2, 0, -1.0),
            (1, 3, 2.0),
            (3, 4, 1.0),
            (4, 5, -1.0),
            (5, 6, 2.0),
            (6, 1, 1.0),
            (3, 7, 1.0),
            (7, 4, 2.0),
        ]);

        let cycle = find_minimum_mean_weight_cycle(&g, ni(0));
        println!("cycle={:?}", cycle);
        assert_eq!(cycle, Some((vec![ni(3), ni(4), ni(5), ni(6), ni(1)], 1.0)));
    }

    ///
    /// graph has a unreachable cycle
    ///
    #[test]
    fn shortest_paths_05() {
        let mut g: DiGraph<(), f64> = DiGraph::new();
        g.extend_with_edges(&[
            // component 1
            (0, 1, 1.0),
            (1, 2, 1.0),
            (0, 2, 1.0),
            // component 2
            (3, 4, 1.0),
            (4, 5, 1.0),
            (5, 6, 1.0),
            (6, 3, 1.0),
        ]);
        // from 0
        let cycle = find_minimum_mean_weight_cycle(&g, ni(0));
        assert_eq!(cycle, None);
        // from 3
        let cycle = find_minimum_mean_weight_cycle(&g, ni(3));
        assert_eq!(cycle, Some((vec![ni(3), ni(4), ni(5), ni(6)], 1.0)));
    }
}
