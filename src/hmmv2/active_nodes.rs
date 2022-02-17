//!
//!
//!
//!
use crate::graph::active_nodes::ActiveNodes;
use crate::hmmv2::common::{PHMMEdge, PHMMModel, PHMMNode};
use crate::vector::graph::NodeVec;
use crate::vector::Storage;
use itertools::Itertools;
use petgraph::graph::NodeIndex;

/// Active nodes methods for PHMMModel and PHMMTable
impl ActiveNodes {
    ///
    /// create a new active nodes consists of our child nodes
    ///
    pub fn to_childs<N: PHMMNode, E: PHMMEdge>(&self, m: &PHMMModel<N, E>) -> ActiveNodes {
        match self {
            ActiveNodes::All => ActiveNodes::All,
            ActiveNodes::Only(nodes) => {
                let childs: Vec<NodeIndex> = nodes
                    .iter()
                    .flat_map(|&node| m.childs(node).map(|(_, child, _)| child))
                    .unique()
                    .collect();
                ActiveNodes::Only(childs)
            }
        }
    }
    ///
    /// create a new active nodes consists of our parent nodes
    ///
    pub fn to_parents<N: PHMMNode, E: PHMMEdge>(&self, m: &PHMMModel<N, E>) -> ActiveNodes {
        match self {
            ActiveNodes::All => ActiveNodes::All,
            ActiveNodes::Only(nodes) => {
                let parents: Vec<NodeIndex> = nodes
                    .iter()
                    .flat_map(|&node| m.parents(node).map(|(_, parent, _)| parent))
                    .unique()
                    .collect();
                ActiveNodes::Only(parents)
            }
        }
    }
}

///
/// NodeVec with ActiveNodes information
///
pub struct ActiveNodeVec<S: Storage> {
    /// NodeVec
    pub vec: NodeVec<S>,
    /// Active Nodes
    pub active_nodes: ActiveNodes,
}

impl<S: Storage> ActiveNodeVec<S> {
    /// Constructor
    pub fn new(
        n_nodes: usize,
        default_value: S::Item,
        active_nodes: ActiveNodes,
    ) -> ActiveNodeVec<S> {
        ActiveNodeVec {
            vec: NodeVec::new(n_nodes, default_value),
            active_nodes,
        }
    }
}
