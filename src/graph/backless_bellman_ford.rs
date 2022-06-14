//!
//! Cusom bellman fords
//!
use super::FloatWeight;
use crate::graph::min_mean_weight_cycle::ShortestPaths;
use petgraph::prelude::*;
use petgraph::visit::{VisitMap, Visitable};

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

//
// tests
//

/*
#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{ei, ni};

    #[test]
    fn shortest_paths_by_edge_01() {
        let mut g: DiGraph<(), f64> = DiGraph::new();
        g.extend_with_edges(&[
            (0, 1, 1.0),
            (1, 2, 4.0),
            (2, 3, 10.0),
            (3, 0, 2.0),
            (3, 4, 2.0),
            (4, 5, 2.0),
            (5, 3, 2.0),
        ]);
        let s = ni(0);
        let p = shortest_paths_by_edge(&g, s, |_, _| true);
        println!("{:?}", p);
        let m = find_minimizer_pair(&g, &p);
        println!("min={:?}", m);
        let (_, e, _) = m.unwrap();
        let path = traceback_preds(&g, e, &p);
        println!("path={:?}", path);
        let cycle = find_cycle(&g, &path);
        println!("cycle={:?}", cycle);
        assert_eq!(cycle, Some(vec![ni(4), ni(5), ni(3)]));

        // prohibit use of e4
        let p = shortest_paths_by_edge(&g, ni(0), |_, t| t != ei(4));
        println!("{:?}", p);
        let m = find_minimizer_pair(&g, &p);
        println!("min={:?}", m);
        let (_, e, _) = m.unwrap();
        let path = traceback_preds(&g, e, &p);
        println!("path={:?}", path);
        let cycle = find_cycle(&g, &path);
        println!("cycle={:?}", cycle);
        assert_eq!(cycle, Some(vec![ni(0), ni(1), ni(2), ni(3)]));
    }
    #[test]
    fn mmwc_edge_01() {
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
        let cycle = find_minimum_mean_weight_cycle(&g, ni(0), |_, _| true);
        println!("{:?}", cycle);
        assert_eq!(cycle, Some((vec![ni(3), ni(4), ni(5), ni(1)], 1.25)));
    }
    #[test]
    fn mmwc_edge_02() {
        let mut g: DiGraph<(), f64> = DiGraph::new();
        g.extend_with_edges(&[
            (0, 1, 3.0),
            (0, 1, 1.0),
            (1, 2, 1.0),
            (2, 0, 1.0),
            (1, 3, 1.0),
            (3, 2, 4.0),
        ]);
        let cycle = find_minimum_mean_weight_cycle(&g, ni(0), |_, _| true);
        println!("{:?}", cycle);
        assert_eq!(cycle, Some((vec![ni(0), ni(1), ni(2)], 1.0)));
    }
    #[test]
    fn mmwc_edge_03() {
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

        let cycle = find_minimum_mean_weight_cycle(&g, ni(0), |_, _| true);
        println!("cycle={:?}", cycle);
        assert_eq!(cycle, Some((vec![ni(2), ni(0), ni(1)], 1.0)));

        let cycle = find_minimum_mean_weight_cycle(&g, ni(0), |_, e| e != ei(1));
        println!("cycle={:?}", cycle);
        assert_eq!(cycle, Some((vec![ni(6), ni(1), ni(3), ni(4), ni(5)], 1.0)));
    }
    #[test]
    fn mmwc_edge_04() {
        let mut g: DiGraph<(), f64> = DiGraph::new();
        g.extend_with_edges(&[
            (0, 1, -5.0),
            (1, 2, -5.0),
            (2, 3, -5.0),
            (3, 4, -5.0),
            (4, 0, -5.0),
            (1, 0, -10.0),
            (2, 1, 10.0),
            (3, 2, 10.0),
            (4, 3, 10.0),
            (0, 4, 10.0),
        ]);

        // (1) mmwc among all cycles
        let cycle = find_minimum_mean_weight_cycle(&g, ni(0), |_, _| true);
        println!("cycle={:?}", cycle);
        assert_eq!(cycle, Some((vec![ni(0), ni(1)], -7.5)));

        // (2) restricted mmwc
        let cycle = find_minimum_mean_weight_cycle(&g, ni(0), |e_a, e_b| {
            e_a.index().abs_diff(e_b.index()) != 5
        });
        println!("cycle={:?}", cycle);
        assert_eq!(cycle, Some((vec![ni(2), ni(3), ni(4), ni(0), ni(1)], -5.0)));
    }
    #[test]
    fn mmwc_edge_05() {
        let mut g: DiGraph<(), f64> = DiGraph::new();
        g.extend_with_edges(&[
            (0, 1, 5.0),
            (1, 2, -5.0),
            (1, 3, 5.0),
            (3, 4, 5.0),
            (4, 5, -100.0),
            (5, 6, 5.0),
            (6, 1, 5.0),
        ]);
        let cycle = find_minimum_mean_weight_cycle(&g, ni(0), |_, _| true);
        println!("cycle={:?}", cycle);
    }
}
*/
