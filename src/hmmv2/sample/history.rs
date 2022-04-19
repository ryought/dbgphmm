//!
//! Struct for storeing sampling results
//!
use super::super::common::{PHMMEdge, PHMMModel, PHMMNode};
use super::super::freq::NodeFreqs;
use super::super::trans_table::EdgeFreqs;
use super::{Emission, State};
use crate::common::{Reads, Sequence};
use crate::graph::genome_graph::GenomeGraphPos;
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
    /// Convert to NodeFreqs (that represents usage frequency of nodes
    /// when emitting the sample)
    ///
    pub fn to_node_freqs<N, E>(&self, phmm: &PHMMModel<N, E>) -> NodeFreqs
    where
        N: PHMMNode,
        E: PHMMEdge,
    {
        let mut nf = NodeFreqs::new(phmm.n_nodes(), 0.0);
        for (state, _) in self.0.iter() {
            match state.to_node_index() {
                Some(v) => nf[v] += 1.0,
                _ => {}
            }
        }
        nf
    }
    /// Convert to EdgeFreqs that represents the usage frequencies of edges
    /// when emitting the sample
    ///
    pub fn to_edge_freqs<N, E>(&self, phmm: &PHMMModel<N, E>) -> EdgeFreqs
    where
        N: PHMMNode,
        E: PHMMEdge,
    {
        let mut ef = EdgeFreqs::new(phmm.n_edges(), 0.0);
        for ((s1, _), (s2, _)) in self.0.iter().tuple_windows() {
            match (s1.to_node_index(), s2.to_node_index()) {
                (Some(v1), Some(v2)) => {
                    // There can be self transition (such as Match(v) -> Ins(v))
                    if v1 != v2 {
                        let e = phmm.find_edge(v1, v2).unwrap();
                        ef[e] += 1.0;
                    }
                }
                _ => {}
            }
        }
        ef
    }
    ///
    /// Convert into genome graph pos list.
    ///
    pub fn to_genome_graph_pos(&self, sg: &SimpleSeqGraph) -> Vec<GenomeGraphPos> {
        self.0
            .iter()
            .filter_map(|(state, _)| {
                state
                    .to_node_index()
                    .map(|v| sg.node_weight(v).unwrap().source())
            })
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

///
/// Vector of multiple historys
///
pub struct Historys(pub Vec<History>);

impl Historys {
    pub fn n_history(&self) -> usize {
        self.0.len()
    }
    pub fn to_sequence(&self, index: usize) -> Sequence {
        self.0[index].to_sequence()
    }
    pub fn to_reads(&self) -> Reads {
        let reads: Vec<Sequence> = self.0.iter().map(|history| history.to_sequence()).collect();
        Reads { reads }
    }
    pub fn iter(&self) -> impl Iterator<Item = &History> + '_ {
        self.0.iter()
    }
}

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
        let s_true = r#"MatchBegin -> -
Match(0) -> A
Ins(1) -> C
Del(2) -> -
Match(3) -> G
End -> -
"#;
        assert_eq!(s, s_true);
        let b = h.to_sequence();
        let b_true = b"ACG";
        assert_eq!(b, b_true);

        assert_eq!(h.total_bases(), 3);
    }

    #[test]
    fn hmm_history_to_freqs() {
        let phmm = mock_linear_phmm(PHMMParams::high_error());
        // choosing seeds..
        for i in 0..20 {
            let hist = phmm.sample(10, i);
            println!("{} {}", i, hist);
        }

        // seed=0 ni(3) -> ni(9)
        let h = phmm.sample(10, 0);
        println!("{}", h);
        let nf = h.to_node_freqs(&phmm);
        assert_abs_diff_eq!(
            nf,
            NodeFreqs::from_slice(&[0.0, 0.0, 0.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0], 0.0)
        );
        let ef = h.to_edge_freqs(&phmm);
        assert_abs_diff_eq!(
            ef,
            EdgeFreqs::from_slice(&[0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], 0.0)
        );

        // seed=2 ni(7) -> ni(9)
        let h = phmm.sample(10, 2);
        println!("{}", h);
        let nf = h.to_node_freqs(&phmm);
        assert_abs_diff_eq!(
            nf,
            NodeFreqs::from_slice(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 0.0)
        );
        let ef = h.to_edge_freqs(&phmm);
        assert_abs_diff_eq!(
            ef,
            EdgeFreqs::from_slice(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0], 0.0)
        );

        // seed=15 ni(0) -> ni(8)
        let h = phmm.sample(10, 15);
        println!("{}", h);
        let nf = h.to_node_freqs(&phmm);
        assert_abs_diff_eq!(
            nf,
            NodeFreqs::from_slice(&[1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0], 0.0)
        );
        let ef = h.to_edge_freqs(&phmm);
        assert_abs_diff_eq!(
            ef,
            EdgeFreqs::from_slice(&[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0], 0.0)
        );
    }
}
