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
    /// Create active_nodes list from NodeVec
    /// by taking the n_active_nodes highest prob nodes from `v`.
    ///
    pub fn from_nodevec<S>(v: &NodeVec<S>, n_active_nodes: usize) -> ActiveNodes
    where
        S: Storage<Item = Prob>,
    {
        let nodes = v
            .iter()
            .sorted_by(|a, b| Ord::cmp(&b.1, &a.1)) // sort by prob in descending order
            .take(n_active_nodes)
            .map(|(node, _)| node)
            .collect();
        ActiveNodes::Only(nodes)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ni;
    use crate::prob::p;
    use crate::vector::dense::DenseStorage;

    #[test]
    fn active_nodes_from_node_vec() {
        let mut v: NodeVec<DenseStorage<Prob>> = NodeVec::new(10, p(0.0));
        v[ni(4)] = p(0.8);
        v[ni(2)] = p(0.6);
        v[ni(1)] = p(0.1);
        v[ni(5)] = p(0.9);
        v[ni(8)] = p(0.5);
        let an = ActiveNodes::from_nodevec(&v, 3);
        assert_eq!(an, ActiveNodes::Only(vec![ni(5), ni(4), ni(2)]));
        println!("{:?}", an);
    }
}
