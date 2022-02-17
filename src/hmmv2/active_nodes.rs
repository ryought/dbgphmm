//!
//!
//!
//!
use super::table::PHMMTable;
use crate::graph::active_nodes::ActiveNodes;
use crate::hmmv2::common::{PHMMEdge, PHMMModel, PHMMNode};
use crate::prob::Prob;
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
    ///
    pub fn fit_to_table<S>(&self, t: &PHMMTable<S>, n_active_nodes: usize) -> ActiveNodes
    where
        S: Storage<Item = Prob>,
    {
        unimplemented!();
    }
    ///
    pub fn fit_to_nodevec<S>(&self, v: &NodeVec<S>, n_active_nodes: usize) -> ActiveNodes
    where
        S: Storage<Item = Prob>,
    {
        unimplemented!();
    }
}
