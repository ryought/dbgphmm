//!
//! Graph constructors
//!
//! From linear seq...
//!
//! From kmer table...
//!

use super::common::{EdgeIndex, NodeIndex, PEdge, PGraph, PModel, PNode};
use crate::common::CopyNum;
use crate::hmm::params::PHMMParams;
use crate::prob::Prob;

pub fn create_linear(seq: &[u8]) -> PModel {
    let mut graph: PGraph = PGraph::new();
    let total_base = seq.len();
    for i in 0..seq.len() {
        graph.add_node(PNode::new(
            1,
            Prob::from_prob(1.0 / total_base as f64),
            true,
            seq[i],
        ));
    }
    for i in 1..seq.len() {
        graph.add_edge(
            NodeIndex::new(i - 1),
            NodeIndex::new(i),
            PEdge::new(Prob::from_prob(1.0)),
        );
    }

    PModel {
        param: PHMMParams::default(),
        graph,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_linear_test() {
        let m = create_linear(b"ATCGGCTAGCT");
        println!("{}", m);
    }
}
