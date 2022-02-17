//!
//! Output related functions of Dbg
//!
use super::dbg::{Dbg, DbgEdge, DbgNode};
use petgraph::dot::Dot;

impl<N, E> std::fmt::Display for Dbg<N, E>
where
    N: DbgNode + std::fmt::Display,
    E: DbgEdge + std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", Dot::with_config(&self.graph, &[]))
    }
}
