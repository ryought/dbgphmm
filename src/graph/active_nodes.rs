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

impl ActiveNodes {
    pub fn is_only(&self) -> bool {
        match self {
            ActiveNodes::Only(_) => true,
            ActiveNodes::All => false,
        }
    }
    pub fn count(&self) -> Option<usize> {
        match self {
            ActiveNodes::Only(nodes) => Some(nodes.len()),
            ActiveNodes::All => None,
        }
    }
}
