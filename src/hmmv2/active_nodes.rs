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
        self.fit_to_nodevec(&t.to_nodevec(), n_active_nodes)
    }
    ///
    /// Create active_nodes list from NodeVec
    /// by taking the n_active_nodes highest prob nodes from `v`.
    ///
    pub fn fit_to_nodevec<S>(&self, v: &NodeVec<S>, n_active_nodes: usize) -> ActiveNodes
    where
        S: Storage<Item = Prob>,
    {
        // v.iter().k_smallest()
        unimplemented!();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ni;
    use crate::prob::p;

    #[test]
    fn itertools_k_smallest() {
        let xs = vec![(p(0.2), ni(0)), (p(0.1), ni(10)), (p(0.5), ni(3))];
        // let sm = xs.iter().k_smallest(2);
        // k_largest
        let ys: Vec<NodeIndex> = xs
            .iter()
            .sorted_by(|a, b| Ord::cmp(&b, &a))
            .take(2)
            .map(|(_, v)| v)
            .copied()
            .collect();
        println!("{:?}", ys);
        assert_eq!(ys, vec![ni(3), ni(0)])
    }
}
