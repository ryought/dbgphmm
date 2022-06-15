//!
//! min_mean_weight with edge conditions
//! custom shortest paths with edge indexing for using edge-adjacency condition
//!
use super::ShortestPaths;
use crate::graph::FloatWeight;
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
    /// Convert edge-indexed `ShortestPathsByEdge` into
    /// node-indexed `ShortestPaths`.
    ///
    /// `nd[k][v]` = (min weight path with k edges from source to v)
    /// `ed[k][e]` = (min weight path with k+1 edges from source ending with edge e)
    ///
    pub fn into_shortest_paths<N, E>(
        self,
        graph: &DiGraph<N, E>,
        source: NodeIndex,
    ) -> ShortestPaths {
        let n = graph.node_count();
        let mut dists = vec![vec![f64::INFINITY; n]; n + 1];
        let mut preds = vec![vec![None; n]; n + 1];

        // assertion of ShortestPathsByEdge
        assert_eq!(self.dists.len(), n);
        assert_eq!(self.preds.len(), n);

        // (1) fill dists_node[0]
        dists[0][source.index()] = 0.0;

        // (2) fill dists_node[k+1] by dists_edge[k]
        for k in 0..n {
            for v in graph.node_indices() {
                // dists[k][v]
                // = min_{e: edge *->v} dists[k-1][e]
                //
                //     e0
                // * -----> v
                let e0 = graph
                    .edges_directed(v, Direction::Incoming)
                    .min_by(|ea, eb| {
                        let da = self.dists[k][ea.id().index()];
                        let db = self.dists[k][eb.id().index()];
                        da.partial_cmp(&db).expect("dists contains nan")
                    });
                match e0 {
                    Some(e0) => {
                        // the min-path (ending with node v) ends with edge e0
                        dists[k + 1][v.index()] = self.dists[k][e0.id().index()];
                        preds[k + 1][v.index()] = Some((e0.source(), e0.id()));
                    }
                    None => {
                        // no incoming edges into the node
                        // so the node is marked as unreachable.
                        dists[k + 1][v.index()] = f64::INFINITY;
                        preds[k + 1][v.index()] = None;
                    }
                }
            }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{ei, ni};

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
        let sp = shortest_paths_by_edge(&g, ni(0), |_, _| true);
        println!("{:?}", sp);
    }
}
