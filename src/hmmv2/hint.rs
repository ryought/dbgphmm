//!
//! Hint information
//!
//! * Mapping
//! * Mappings
//!
use super::common::{PHMMEdge, PHMMModel, PHMMNode};
use super::table::{PHMMOutput, MAX_ACTIVE_NODES};
use crate::common::{ReadCollection, Seq};
use crate::prob::Prob;
use crate::utils::progress_common_style;
use arrayvec::ArrayVec;
use fnv::FnvHashMap as HashMap;
use indicatif::ParallelProgressIterator;
use itertools::{izip, Itertools};
use petgraph::graph::NodeIndex;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// Mapping: Hint/Cache of mapping for a emission sequence
///
/// Define candidate nodes for each emission
///
/// `emission[i]` is emitted by nodes in `self.nodes[i]` with probabilities `self.probs[i]`
///
#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
pub struct Mapping {
    pub nodes: Vec<ArrayVec<NodeIndex, MAX_ACTIVE_NODES>>,
    pub probs: Vec<ArrayVec<Prob, MAX_ACTIVE_NODES>>,
}

impl Mapping {
    ///
    /// Constructor from vec of vecs
    ///
    pub fn from_nodes_and_probs(vs: &Vec<ArrayVec<(NodeIndex, Prob), MAX_ACTIVE_NODES>>) -> Self {
        Mapping {
            nodes: vs
                .iter()
                .map(|v| {
                    let mut vec = ArrayVec::new();
                    for (node, _) in v.iter() {
                        vec.push(*node);
                    }
                    vec
                })
                .collect(),
            probs: vs
                .iter()
                .map(|v| {
                    let mut vec = ArrayVec::new();
                    for (_, prob) in v.iter() {
                        vec.push(*prob);
                    }
                    vec
                })
                .collect(),
        }
    }
    ///
    /// Given `node_map: node -> corresponding nodes`
    ///
    /// # TODO
    /// * depends on n_active_nodes
    ///
    pub fn map_nodes<'a, F>(mut self, node_map: F) -> Self
    where
        F: Fn(NodeIndex) -> ArrayVec<NodeIndex, MAX_ACTIVE_NODES>,
    {
        let mut ret = Vec::new();
        for (ns, ps) in izip!(&self.nodes, &self.probs) {
            // let mut m: HashMap<NodeIndex, Prob> = HashMap::default();
            // for (&node, &prob) in izip!(ns, ps) {
            //     let nodes_after = node_map(node);
            //     for &node_after in nodes_after.iter() {
            //         *m.entry(node_after).or_insert(Prob::zero()) += prob / nodes_after.len();
            //     }
            // }

            let vs: ArrayVec<_, MAX_ACTIVE_NODES> = izip!(ns, ps)
                .into_iter()
                .flat_map(|(&node, &prob)| {
                    let nodes_after = node_map(node);
                    let c = nodes_after.len();
                    nodes_after
                        .into_iter()
                        .map(|node_after| (node_after, prob))
                        .collect::<Vec<_>>()
                })
                // .sorted_by_key(|(_, prob)| *prob)
                // .rev()
                .take(MAX_ACTIVE_NODES)
                .collect();
            ret.push(vs);
        }
        Mapping::from_nodes_and_probs(&ret)
    }
    ///
    /// Get candidate nodes of `emissions[index]`
    ///
    pub fn nodes(&self, index: usize) -> &[NodeIndex] {
        &self.nodes[index]
    }
    ///
    /// Length of the read (emission sequence)
    ///
    pub fn len(&self) -> usize {
        self.nodes.len()
    }
    ///
    /// String representations of nodes and probs of `emissions[index]`
    ///
    pub fn to_nodes_string(&self, index: usize) -> String {
        izip!(&self.nodes[index], &self.probs[index])
            .map(|(n, p)| format!("{}:{:.1}", n.index(), p.to_value() * 100.0))
            .join(",")
    }
}

impl std::fmt::Display for Mapping {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for i in 0..self.len() {
            writeln!(f, "{}\t{}", i, self.to_nodes_string(i))?;
        }
        Ok(())
    }
}

impl PHMMOutput {
    ///
    /// Create Mapping of this emission sequence
    ///
    pub fn to_mapping(&self, n_active_nodes: usize) -> Mapping {
        let ret = self
            .iter_emit_probs()
            .skip(1)
            .map(|state_probs| state_probs.top_nodes_with_prob(n_active_nodes))
            .collect();
        Mapping::from_nodes_and_probs(&ret)
    }
}

///
/// Mapping information for reads
///
#[derive(Clone, Debug, PartialEq)]
pub struct Mappings(Vec<Mapping>);

impl Mappings {
    pub fn new(v: Vec<Mapping>) -> Self {
        Mappings(v)
    }
}

impl std::ops::Index<usize> for Mappings {
    type Output = Mapping;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl IntoIterator for Mappings {
    type Item = Mapping;
    type IntoIter = std::vec::IntoIter<Mapping>;
    fn into_iter(self) -> std::vec::IntoIter<Mapping> {
        self.0.into_iter()
    }
}

impl<'a> IntoIterator for &'a Mappings {
    type Item = &'a Mapping;
    type IntoIter = std::slice::Iter<'a, Mapping>;
    fn into_iter(self) -> std::slice::Iter<'a, Mapping> {
        self.0.iter()
    }
}

impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    /// Generate Mappings = Vec<Mapping> (Mapping for all reads) using the model.
    ///
    pub fn generate_mappings<S: Seq>(
        &self,
        reads: &ReadCollection<S>,
        mappings: Option<&Mappings>,
    ) -> Mappings {
        Mappings(
            reads
                .par_iter()
                .enumerate()
                .progress_with_style(progress_common_style())
                .map(|(i, seq)| {
                    let output = if let Some(mappings) = mappings {
                        self.run_with_mapping(seq.as_ref(), &mappings[i])
                    } else {
                        self.run_sparse(seq.as_ref())
                    };
                    output.to_mapping(self.param.n_active_nodes)
                })
                .collect(),
        )
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ni;
    use crate::prob::p;

    #[test]
    fn mapping_node_convert() {
        let m = Mapping {
            nodes: vec![
                vec![ni(0), ni(1)], // emission[0]
                vec![ni(2), ni(3)], // emission[1]
            ],
            probs: vec![
                vec![p(0.6), p(0.4)], // emission[0]
                vec![p(0.9), p(0.1)], // emission[0]
            ],
        };
        println!("before");
        println!("{}", m);

        // case1
        // v -> [v+1]
        let m_1 = m.map_nodes(|node| vec![ni(node.index() + 1)]);
        println!("after");
        println!("{}", m_1);
        let m_1_true = Mapping {
            nodes: vec![
                vec![ni(1), ni(2)], // emission[0]
                vec![ni(3), ni(4)], // emission[1]
            ],
            probs: vec![
                vec![p(0.6), p(0.4)], // emission[0]
                vec![p(0.9), p(0.1)], // emission[0]
            ],
        };
        assert_eq!(m_1, m_1_true);

        // case2
        // v -> [v, v+1]
        let m_2 = m.map_nodes(|node| vec![node, ni(node.index() + 1)]);
        println!("after");
        println!("{}", m_2);
    }
}
