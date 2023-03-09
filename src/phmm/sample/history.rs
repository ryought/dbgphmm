//!
//! Struct for storeing sampling results
//!
use super::super::{PHMMEdge, PHMMNode, PHMM};
use super::{Emission, State};
use crate::common::{PositionedReads, PositionedSequence, Reads, Sequence};
use crate::graph::genome_graph::{GenomeGraphPos, GenomeGraphPosVec};
use crate::graph::seq_graph::SimpleSeqGraph;
use itertools::Itertools;

///
/// Struct for storing sampling results from HMM
///
pub struct History(pub Vec<(State, Emission)>);

impl History {
    ///
    /// Constructor of empty sample store.
    ///
    pub fn new() -> Self {
        History(Vec::new())
    }
    ///
    /// Total history length
    /// This can be different from the length of emitted bases
    ///
    pub fn len(&self) -> usize {
        self.0.len()
    }
    ///
    /// Append a new state and its emission
    ///
    pub fn push(&mut self, state: State, emission: Emission) {
        self.0.push((state, emission));
    }
    ///
    /// Get a total bases count of the emissions
    ///
    pub fn total_bases(&self) -> usize {
        self.0
            .iter()
            .filter(|(_, emission)| emission.is_base())
            .count()
    }
    ///
    /// Create base sequence `Vec<u8>` from sampling history
    ///
    pub fn to_sequence(&self) -> Sequence {
        self.0
            .iter()
            .filter_map(|(_, emission)| match emission {
                Emission::Base(base) => Some(base),
                Emission::Empty => None,
            })
            .copied()
            .collect()
    }
}

//
// Display
//
impl std::fmt::Display for History {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        for (state, emission) in self.0.iter() {
            writeln!(f, "{} -> {}", state, emission)?;
        }
        Ok(())
    }
}

// use super::super::freq::NodeFreqs;
// use super::super::trans_table::EdgeFreqs;
// impl History {
//     /// Convert to NodeFreqs (that represents usage frequency of nodes
//     /// when emitting the sample)
//     ///
//     pub fn to_node_freqs<N, E>(&self, phmm: &PHMMModel<N, E>) -> NodeFreqs
//     where
//         N: PHMMNode,
//         E: PHMMEdge,
//     {
//         let mut nf = NodeFreqs::new(phmm.n_nodes(), 0.0);
//         for (state, _) in self.0.iter() {
//             match state.to_node_index() {
//                 Some(v) => nf[v] += 1.0,
//                 _ => {}
//             }
//         }
//         nf
//     }
//     /// Convert to EdgeFreqs that represents the usage frequencies of edges
//     /// when emitting the sample
//     ///
//     pub fn to_edge_freqs<N, E>(&self, phmm: &PHMMModel<N, E>) -> EdgeFreqs
//     where
//         N: PHMMNode,
//         E: PHMMEdge,
//     {
//         let mut ef = EdgeFreqs::new(phmm.n_edges(), 0.0);
//         for ((s1, _), (s2, _)) in self.0.iter().tuple_windows() {
//             match (s1.to_node_index(), s2.to_node_index()) {
//                 (Some(v1), Some(v2)) => {
//                     // There can be self transition (such as Match(v) -> Ins(v))
//                     if v1 != v2 {
//                         let e = phmm.find_edge(v1, v2).unwrap();
//                         ef[e] += 1.0;
//                     }
//                 }
//                 _ => {}
//             }
//         }
//         ef
//     }
// }

//
// tests
//
#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ni;
    use crate::hmmv2::mocks::mock_linear_phmm;
    use crate::hmmv2::params::PHMMParams;

    #[test]
    fn hmm_history_to_sequence() {
        let mut h = History::new();
        h.push(State::MatchBegin, Emission::Empty);
        h.push(State::Match(ni(0)), Emission::Base(b'A'));
        h.push(State::Ins(ni(1)), Emission::Base(b'C'));
        h.push(State::Del(ni(2)), Emission::Empty);
        h.push(State::Match(ni(3)), Emission::Base(b'G'));
        h.push(State::End, Emission::Empty);
        let s = h.to_string();
        println!("{}", s);
        let s_true = "MB -> -\nM(n0) -> A\nI(n1) -> C\nD(n2) -> -\nM(n3) -> G\nE -> -\n";
        assert_eq!(s, s_true);
        let b = h.to_sequence();
        let b_true = b"ACG";
        assert_eq!(b, b_true);

        assert_eq!(h.total_bases(), 3);
    }
}
