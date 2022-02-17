use super::{EDbg, EDbgEdge, EDbgNode};
use petgraph::dot::Dot;

///
/// Normal DOT output
///
impl<N, E> std::fmt::Display for EDbg<N, E>
where
    N: EDbgNode + std::fmt::Display,
    E: EDbgEdge + std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", Dot::with_config(&self.graph, &[]))
    }
}

//
// TODO cytoscape output
//
