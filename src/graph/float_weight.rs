//!
//! FloatWeight
//!
use petgraph::prelude::*;

pub trait FloatWeight {
    fn float_weight(&self) -> f64;
    // fn epsilon(&self) -> f64;
}

impl FloatWeight for f64 {
    fn float_weight(&self) -> f64 {
        *self
    }
}

///
/// Calculate total weight of path (a list of edges)
///
pub fn total_weight<N, E: FloatWeight>(graph: &DiGraph<N, E>, edges: &[EdgeIndex]) -> f64 {
    edges
        .iter()
        .map(|&e| {
            let ew = graph.edge_weight(e).unwrap();
            ew.float_weight()
        })
        .sum()
}

///
/// Determine if a cycle given by edges is a negative cycle or not.
///
pub fn is_negative_cycle<N, E: FloatWeight>(graph: &DiGraph<N, E>, edges: &[EdgeIndex]) -> bool {
    total_weight(graph, edges) < 0.0
}

/// Find the minimum weight edge among all parallel edges between v and w
/// Input: two nodes (v,w) in a graph
/// Output: minimum weight edge among all parallel edge (v,w)
///
/// (Used in `node_list_to_edge_list`)
fn pick_minimum_weight_edge<N, E: FloatWeight>(
    graph: &DiGraph<N, E>,
    v: NodeIndex,
    w: NodeIndex,
) -> EdgeIndex {
    let er = graph
        .edges_connecting(v, w)
        .min_by(|e1, e2| {
            let w1 = e1.weight().float_weight();
            let w2 = e2.weight().float_weight();
            w1.partial_cmp(&w2).unwrap()
        })
        .unwrap();
    let e = er.id();
    e
}

///
/// Convert "a cycle as nodes [NodeIndex]" into "a cycle as edges [EdgeIndex]",
/// by choosing the minimum weight edge if there are parallel edges
///
pub fn node_list_to_edge_list<N, E: FloatWeight>(
    graph: &DiGraph<N, E>,
    nodes: &[NodeIndex],
) -> Vec<EdgeIndex> {
    let mut edges = Vec::new();
    let n = nodes.len();

    // convert (nodes[i], nodes[i+1]) into an edge
    for i in 0..n {
        let v = nodes[i];
        let w = nodes[(i + 1) % n];
        let edge = pick_minimum_weight_edge(graph, v, w);
        edges.push(edge);
    }

    edges
}

///
/// Convert edge list into node list
///
pub fn edge_list_to_node_list<N, E>(graph: &DiGraph<N, E>, edges: &[EdgeIndex]) -> Vec<NodeIndex> {
    let mut nodes = Vec::new();

    // (1) add first edge
    let (v0, _) = graph
        .edge_endpoints(edges[0])
        .expect("edge is not in the graph");
    nodes.push(v0);

    // (2) add target node of each edge
    for &edge in edges.iter() {
        let (_, v) = graph
            .edge_endpoints(edge)
            .expect("edge is not in the graph");
        nodes.push(v);
    }

    nodes
}

///
/// Convert cycle of edge into node list
///
pub fn edge_cycle_to_node_cycle<N, E>(
    graph: &DiGraph<N, E>,
    cycle: &[EdgeIndex],
) -> Vec<NodeIndex> {
    cycle
        .iter()
        .map(|&edge| {
            let (v, _) = graph
                .edge_endpoints(edge)
                .expect("edge is not in the graph");
            v
        })
        .collect()
}
