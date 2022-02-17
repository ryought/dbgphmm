//!
//! ActiveNodes store
//!
use petgraph::graph::NodeIndex;

///
/// Collection of active nodes
///
#[derive(Debug, Clone, PartialEq)]
pub enum ActiveNodes {
    /// All nodes are active
    All,
    /// Only nodes specified are active
    Only(Vec<NodeIndex>),
}
