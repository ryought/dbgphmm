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
