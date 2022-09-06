//!
//! Cycle in graph
//!
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

///
/// cycle (as a list of edges)
///
#[derive(Debug, Clone)]
pub struct Cycle(Vec<EdgeIndex>);

impl Cycle {
    /// constructor from vec of edgeindex
    pub fn new(edges: Vec<EdgeIndex>) -> Cycle {
        Cycle(edges)
    }
    /// normalize cycle
    pub fn normalize(self) -> Cycle {
        unimplemented!();
    }
}
