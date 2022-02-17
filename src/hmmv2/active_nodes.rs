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
use itertools::{chain, Itertools};
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
    /// create a new active nodes consists of our parent nodes and current nodes
    ///
    pub fn to_parents_and_us<N: PHMMNode, E: PHMMEdge>(&self, m: &PHMMModel<N, E>) -> ActiveNodes {
        match self {
            ActiveNodes::All => ActiveNodes::All,
            ActiveNodes::Only(nodes) => {
                let parents: Vec<NodeIndex> = nodes
                    .iter()
                    .flat_map(|&node| m.parents(node).map(|(_, parent, _)| parent))
                    .chain(nodes.iter().copied())
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
    ///
    /// merge another active nodes into self.
    ///
    /// * if self or other is All, merged nodes is also All
    /// * if both self and other are Only, nodes will be merged.
    ///
    pub fn merge(&self, other: &ActiveNodes) -> ActiveNodes {
        match (self, other) {
            (ActiveNodes::Only(n1), ActiveNodes::Only(n2)) => {
                let n = chain!(n1.iter().copied(), n2.iter().copied())
                    .unique()
                    .collect();
                ActiveNodes::Only(n)
            }
            _ => ActiveNodes::All,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ni;
    use crate::hmmv2::mocks::mock_linear_phmm;
    use crate::hmmv2::params::PHMMParams;
    use crate::prob::p;
    use crate::vector::dense::DenseStorage;

    #[test]
    fn active_nodes_update() {
        let phmm = mock_linear_phmm(PHMMParams::zero_error());
        let an = ActiveNodes::Only(vec![ni(1)]);
        assert_eq!(an.to_childs(&phmm), ActiveNodes::Only(vec![ni(2)]));
        assert_eq!(an.to_parents(&phmm), ActiveNodes::Only(vec![ni(0)]));
        assert_eq!(
            an.to_parents_and_us(&phmm),
            ActiveNodes::Only(vec![ni(0), ni(1)])
        );

        let an = ActiveNodes::Only(vec![ni(0)]);
        assert_eq!(an.to_parents(&phmm), ActiveNodes::Only(vec![]));

        let an = ActiveNodes::Only(vec![ni(3), ni(6), ni(6)]);
        assert_eq!(an.to_parents(&phmm), ActiveNodes::Only(vec![ni(2), ni(5)]));
    }

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
    #[test]
    fn active_nodes_merge() {
        let an1 = ActiveNodes::Only(vec![ni(1), ni(2)]);
        let an2 = ActiveNodes::Only(vec![ni(1), ni(3), ni(5)]);
        let an = an1.merge(&an2);
        println!("{:?}", an1);
        println!("{:?}", an2);
        println!("{:?}", an);
        assert_eq!(an, ActiveNodes::Only(vec![ni(1), ni(2), ni(3), ni(5)]));
    }
}
