use petgraph::graph::{DiGraph, EdgeIndex, Graph, NodeIndex, UnGraph};

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
