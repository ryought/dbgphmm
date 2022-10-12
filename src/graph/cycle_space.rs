//!
//! cycle space
//!
//! enumerate all cycles in the undirected graph
//!
use super::cycle::{Cycle, SimpleCycle};
use petgraph::unionfind::UnionFind;

#[derive(Clone, Debug)]
struct CycleSpace {
    ///
    /// a set of all cycle basis.
    ///
    basis: Vec<SimpleCycle>,
}
