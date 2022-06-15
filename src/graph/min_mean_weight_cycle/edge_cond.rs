//!
//! min_mean_weight with edge conditions
//! custom shortest paths with edge indexing for using edge-adjacency condition
//!
use super::ShortestPaths;
use crate::graph::float_weight::total_weight;
use crate::graph::FloatWeight;
use fnv::FnvHashMap as HashMap;
use itertools::Itertools;
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

        // TODO fix assertion of ShortestPathsByEdge
        // assert_eq!(self.dists.len(), n);
        // assert_eq!(self.preds.len(), n);

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
    let n = graph.edge_count();
    let ix = |edge: EdgeIndex| edge.index();
    let eps = f64::EPSILON;

    // (1) Initialize
    //     e
    // s ----> *
    let mut dists = vec![vec![f64::INFINITY; n]; n + 1];
    let mut preds = vec![vec![None; n]; n + 1];
    for edge in graph.edges_directed(source, Direction::Outgoing) {
        dists[0][ix(edge.id())] = edge.weight().float_weight();
    }

    // (2) Update
    for k in 1..=n {
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
                        preds[k][ix(e)] = Some(ep);
                    }
                }
            }
        }
    }

    ShortestPathsByEdge { dists, preds }
}

///
/// Find a minimizer pair `(k, e)`
/// that satisfies `min_e max_k (F[n][e] - F[k][e]) / (n-k)`.
///
fn find_minimizer_pair<N, E: FloatWeight>(
    graph: &DiGraph<N, E>,
    paths: &ShortestPathsByEdge,
) -> Option<(usize, EdgeIndex, f64)> {
    let n = graph.edge_count();
    (0..n)
        .filter_map(|e| {
            (0..n)
                .filter_map(|k| {
                    let fnv = paths.dists[n][e];
                    let fkv = paths.dists[k][e];
                    if fnv != f64::INFINITY && fkv != f64::INFINITY {
                        Some((k, (fnv - fkv) / (n - k) as f64))
                    } else {
                        None
                    }
                })
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                .map(|(k, score)| (k, EdgeIndex::new(e), score))
        })
        .min_by(|(_, _, a), (_, _, b)| a.partial_cmp(b).unwrap())
}

///
/// Traceback a path from source to target
///
pub fn traceback<N, E: FloatWeight>(
    graph: &DiGraph<N, E>,
    target: EdgeIndex,
    paths: &ShortestPathsByEdge,
) -> Vec<EdgeIndex> {
    let n = graph.edge_count();
    let ix = |edge: EdgeIndex| edge.index();

    let mut edge = target;
    let mut edges = vec![target];

    for k in (1..=n).rev() {
        let edge_pred = match paths.preds[k][ix(edge)] {
            Some(e) => e,
            None => panic!("no parent"),
        };
        edges.push(edge_pred);
        edge = edge_pred;
    }

    edges.reverse();
    edges
}

///
///
pub fn find_all_cycles_in_path<N, E>(
    _graph: &DiGraph<N, E>,
    path: &[EdgeIndex],
) -> Vec<(usize, usize)> {
    // create occurrence table
    let mut pos: HashMap<EdgeIndex, Vec<usize>> = HashMap::default();
    for (i, &edge) in path.iter().enumerate() {
        pos.entry(edge).or_insert_with(|| Vec::new()).push(i);
    }

    // for each edge that is used multiple times
    let mut cycles = Vec::new();
    for (_edge, occ) in &pos {
        if occ.len() > 1 {
            for (&i, &j) in occ.iter().sorted().tuple_combinations() {
                cycles.push((i, j));
            }
        }
    }

    cycles
}

///
/// Find a cycle in a path
///
pub fn find_min_mean_weight_cycle_in_path<N, E: FloatWeight>(
    graph: &DiGraph<N, E>,
    path: &[EdgeIndex],
) -> Option<Vec<EdgeIndex>> {
    let cycles = find_all_cycles_in_path(graph, path);

    cycles
        .iter()
        .min_by(|(ia, ja), (ib, jb)| {
            let wa = total_weight(graph, &path[*ia..*ja]) / (ja - ia) as f64;
            let wb = total_weight(graph, &path[*ib..*jb]) / (jb - ib) as f64;
            wa.partial_cmp(&wb).expect("dists contains nan")
        })
        .map(|(i, j)| path[*i..*j].to_vec())
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{ei, ni};

    fn into_path(edges: &[usize]) -> Vec<EdgeIndex> {
        edges.iter().map(|&i| ei(i)).collect()
    }

    #[test]
    fn find_all_cycles_00() {
        let mut g: DiGraph<(), f64> = DiGraph::new();
        g.extend_with_edges(&[
            // cycle 1
            (0, 1, 1.0),
            (1, 2, 1.0),
            (2, 3, 1.0),
            (3, 0, 1.0),
            // cycle 2
            (2, 4, 1.0),
            (4, 5, 1.0),
            (5, 6, 1.0),
            (6, 2, 1.0),
        ]);
        // (a)
        let ix = vec![0, 1, 2, 3, 0];
        let cycles = find_all_cycles_in_path(&g, &into_path(&ix));
        println!("{:?}", cycles);
        let cycle = find_min_mean_weight_cycle_in_path(&g, &into_path(&ix));
        println!("{:?}", cycle);
        assert_eq!(cycle, Some(into_path(&[0, 1, 2, 3])));

        // (b)
        let ix = vec![0, 1, 4, 5, 6, 7, 4, 5, 6, 7, 2, 3, 0];
        let cycles = find_all_cycles_in_path(&g, &into_path(&ix));
        println!("{:?}", cycles);
        let cycle = find_min_mean_weight_cycle_in_path(&g, &into_path(&ix));
        println!("{:?}", cycle);
        assert_eq!(cycle, Some(into_path(&[5, 6, 7, 4])));
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
        let sp = shortest_paths_by_edge(&g, ni(0), |_, _| true);
        println!("{:?}", sp);
        let mp = find_minimizer_pair(&g, &sp);
        println!("{:?}", mp);
    }
}
