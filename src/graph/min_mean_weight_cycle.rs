//!
//! Find a minimum mean-weight cycle in petgraph::DiGraph
//!
//! Karp's minimum mean-weight cycle algorithm
//!
//! # References
//!
//! * https://walkccc.me/CLRS/Chap24/Problems/24-5/
//!
use petgraph::prelude::*;

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
    /// `preds[k: length of path][v: node] = w: node`
    /// means that "min weight path `F[k][v]` ends with the `w->v` edge"
    ///
    preds: Vec<Vec<Option<NodeIndex>>>,
}

impl ShortestPaths {}

pub trait FloatWeight {
    fn float_weight(&self) -> f64;
}

impl FloatWeight for f64 {
    fn float_weight(&self) -> f64 {
        *self
    }
}

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
                // an edge u->v with weight w
                let v = edge.target();
                let w = edge.weight().float_weight();
                // if s->u->v is shorter than s->v, update the route.
                if dists[k - 1][ix(u)] + w < dists[k][ix(v)] {
                    dists[k][ix(v)] = dists[k - 1][ix(u)] + w;
                    preds[k][ix(v)] = Some(u);
                }
            }
        }
    }

    ShortestPaths { dists, preds }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ni;

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
        assert_eq!(sp.preds[1][1], Some(ni(0)));
        // 0->2
        assert_eq!(sp.dists[1][2], 1.0);
        assert_eq!(sp.preds[1][2], Some(ni(0)));
        // 0->2->3
        assert_eq!(sp.dists[2][3], 1.0 + 1.0);
        assert_eq!(sp.preds[2][3], Some(ni(2)));
        // 0->2->4
        assert_eq!(sp.dists[2][4], 1.0 + 2.0);
        assert_eq!(sp.preds[2][4], Some(ni(2)));
    }
}
