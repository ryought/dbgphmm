//!
//! Hint information
//!
//! * Mapping
//! * Mappings
//!
use super::common::{PHMMEdge, PHMMModel, PHMMNode};
use super::freq::NodeFreqs;
use super::table::{PHMMOutput, MAX_ACTIVE_NODES};
use crate::common::{ReadCollection, Seq};
use crate::prob::Prob;
use crate::utils::progress_common_style;
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
    pub nodes: Vec<Vec<NodeIndex>>,
    pub probs: Vec<Vec<Prob>>,
}

impl Mapping {
    ///
    /// Constructor from vec of vecs
    ///
    pub fn from_nodes_and_probs(vs: &Vec<Vec<(NodeIndex, Prob)>>) -> Self {
        Mapping {
            nodes: vs
                .iter()
                .map(|v| {
                    let mut vec = Vec::new();
                    for (node, _) in v.iter() {
                        vec.push(*node);
                    }
                    vec
                })
                .collect(),
            probs: vs
                .iter()
                .map(|v| {
                    let mut vec = Vec::new();
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
    pub fn map_nodes<'a, F>(&self, node_map: F) -> Self
    where
        F: Fn(NodeIndex) -> Vec<NodeIndex>,
    {
        let mut ret = Vec::new();
        for (ns, ps) in izip!(&self.nodes, &self.probs) {
            let mut m: HashMap<NodeIndex, Prob> = HashMap::default();
            for (&node, &prob) in izip!(ns, ps) {
                let nodes_after = node_map(node);
                for &node_after in nodes_after.iter() {
                    *m.entry(node_after).or_insert(Prob::zero()) += prob / nodes_after.len();
                }
            }
            let vs: Vec<_> = m
                .into_iter()
                .sorted_by_key(|(_, prob)| *prob)
                .rev()
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
    ///
    /// Create Mapping of this emission sequence
    ///
    pub fn to_mapping_by_score_ratio(&self, max_ratio: f64) -> Mapping {
        let ret = self
            .iter_emit_probs()
            .skip(1)
            .map(|state_probs| state_probs.top_nodes_with_prob_by_score_ratio(max_ratio))
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
    pub fn n_reads(&self) -> usize {
        self.0.len()
    }
    ///
    /// Convert Mappings into NodeFreqs, the expected usage of each node.
    ///
    pub fn to_node_freqs(&self, n_nodes: usize) -> NodeFreqs {
        let mut freqs: NodeFreqs = NodeFreqs::new(n_nodes, 0.0, true);
        for mapping in self.into_iter() {
            for (ns, ps) in izip!(&mapping.nodes, &mapping.probs) {
                for (&n, &p) in izip!(ns, ps) {
                    freqs[n] += p.to_value();
                }
            }
        }
        freqs
    }
}

impl std::ops::Index<usize> for Mappings {
    type Output = Mapping;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
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
        use_max_ratio: bool,
    ) -> Mappings {
        Mappings(
            reads
                .par_iter()
                // .iter()
                .enumerate()
                .progress_with_style(progress_common_style())
                .map(|(i, seq)| {
                    let output = if let Some(mappings) = mappings {
                        self.run_with_mapping(seq.as_ref(), &mappings[i])
                    } else {
                        self.run_sparse_adaptive(seq.as_ref(), use_max_ratio)
                    };

                    if use_max_ratio {
                        output.to_mapping_by_score_ratio(self.param.active_node_max_ratio)
                    } else {
                        output.to_mapping(self.param.n_active_nodes)
                    }
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
