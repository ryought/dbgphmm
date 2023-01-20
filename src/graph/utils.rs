use fnv::FnvHashMap as HashMap;
use petgraph::graph::{DiGraph, EdgeIndex, Graph, NodeIndex, UnGraph};
use petgraph::Direction;

pub fn to_node_list<N: Clone, E>(graph: &DiGraph<N, E>) -> Vec<(usize, N)> {
    graph
        .node_indices()
        .map(|node| {
            let weight = graph.node_weight(node).unwrap().clone();
            (node.index(), weight)
        })
        .collect()
}

pub fn to_edge_list<N, E: Clone>(graph: &DiGraph<N, E>) -> Vec<(usize, usize, usize, E)> {
    graph
        .edge_indices()
        .map(|edge| {
            let weight = graph.edge_weight(edge).unwrap().clone();
            let (source, target) = graph.edge_endpoints(edge).unwrap();
            (edge.index(), source.index(), target.index(), weight)
        })
        .collect()
}

///
/// count the number of (in_degree, out_degree) of nodes.
///
pub fn degree_stats<N, E>(graph: &DiGraph<N, E>) -> HashMap<(usize, usize), usize> {
    let mut h = HashMap::default();
    for node in graph.node_indices() {
        let in_degree = graph.edges_directed(node, Direction::Incoming).count();
        let out_degree = graph.edges_directed(node, Direction::Outgoing).count();
        *h.entry((in_degree, out_degree)).or_insert(0) += 1;
    }
    h
}
