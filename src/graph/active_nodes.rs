use petgraph::graph::NodeIndex;
///
/// Collection of active nodes
///
#[derive(Debug, Clone)]
pub enum ActiveNodes {
  /// All nodes are active
  All,
  /// Only nodes specified are active
  Only(Vec<NodeIndex>),
}
