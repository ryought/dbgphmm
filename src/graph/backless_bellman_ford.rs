//!
//! Cusom bellman fords
//!
use super::FloatWeight;
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
    let mut dists = vec![vec![f64::INFINITY; n_edges]; n_nodes + 1];
    let mut preds = vec![vec![None; n_edges]; n_nodes + 1];
    for edge in graph.edges_directed(source, Direction::Outgoing) {
        dists[0][ix(edge.id())] = edge.weight().float_weight();
    }

    // (2) Update
    for k in 1..=n_nodes {
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

fn find_minimizer_pair<N, E: FloatWeight>(
    graph: &DiGraph<N, E>,
    paths: &ShortestPathsByEdge,
) -> Option<(usize, EdgeIndex, f64)> {
    let n = graph.node_count();
    let ix = |edge: EdgeIndex| edge.index();
    graph
        .edge_indices()
        .filter_map(|e| {
            (0..n)
                .filter_map(|k| {
                    let fnv = paths.dists[n][ix(e)];
                    let fkv = paths.dists[k][ix(e)];
                    if fnv != f64::INFINITY && fkv != f64::INFINITY {
                        Some((k, (fnv - fkv) / (n - k) as f64))
                    } else {
                        None
                    }
                })
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                .map(|(k, score)| (k, e, score))
        })
        .min_by(|(_, _, a), (_, _, b)| a.partial_cmp(b).unwrap())
}

///
/// Traceback a path from source to target
/// by using `preds` in paths
///
fn traceback_preds<N, E: FloatWeight>(
    graph: &DiGraph<N, E>,
    edge_last: EdgeIndex,
    paths: &ShortestPathsByEdge,
) -> Vec<EdgeIndex> {
    let n = graph.node_count();
    let ix = |edge: EdgeIndex| edge.index();

    let mut edge = edge_last;
    let mut path = vec![edge_last];

    for k in (1..=n).rev() {
        let pred = match paths.preds[k][ix(edge)] {
            Some(pred) => pred,
            None => panic!("no parent"),
        };
        path.push(pred);
        edge = pred;
    }

    path.reverse();
    path
}

///
/// Find a cycle in a path (EdgeList)
///
fn find_cycle<N, E>(graph: &DiGraph<N, E>, path: &[EdgeIndex]) -> Option<Vec<NodeIndex>> {
    let mut v = graph.visit_map();
    let mut cycle = Vec::new();

    for &edge in path.iter().rev() {
        let (_, node) = graph
            .edge_endpoints(edge)
            .expect("an edge in the given path does not exist in the graph");
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
pub fn find_minimum_mean_weight_cycle<N, E, F>(
    graph: &DiGraph<N, E>,
    source: NodeIndex,
    edge_moveable: F,
) -> Option<(Vec<NodeIndex>, f64)>
where
    E: FloatWeight,
    F: Fn(EdgeIndex, EdgeIndex) -> bool,
{
    let sp = shortest_paths_by_edge(graph, source, edge_moveable);
    match find_minimizer_pair(graph, &sp) {
        Some((_, v, mean_weight)) => {
            let path = traceback_preds(graph, v, &sp);
            let cycle = find_cycle(graph, &path)
                .expect("minimizer pair was found, but no cycle was found when tracebacking");
            Some((cycle, mean_weight))
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
}
