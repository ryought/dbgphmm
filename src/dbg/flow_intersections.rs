use super::dbg::{Dbg, DbgEdge, DbgNode};
use crate::dbg::edge_centric::SimpleEDbg;
use crate::em::extension::flow_intersection::FlowIntersection;
use crate::hmmv2::trans_table::EdgeFreqs;
use crate::kmer::kmer::KmerLike;
use petgraph::graph::NodeIndex;

pub struct FlowIntersections<'a, K: KmerLike> {
    edbg: SimpleEDbg<K>,
    current_node: usize,
    freqs: &'a EdgeFreqs,
}

impl<'a, K: KmerLike> Iterator for FlowIntersections<'a, K> {
    type Item = FlowIntersection<K>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.current_node < self.edbg.n_nodes() {
            let intersection = self
                .edbg
                .flow_intersection(NodeIndex::new(self.current_node), self.freqs);
            self.current_node += 1;
            Some(intersection)
        } else {
            None
        }
    }
}

///
/// Intersections
///
impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    ///
    /// get an iterator over intersections (k-1-mer)
    ///
    pub fn iter_flow_intersections<'a>(
        &'a self,
        freqs: &'a EdgeFreqs,
    ) -> FlowIntersections<'a, N::Kmer> {
        FlowIntersections {
            edbg: self.to_edbg(),
            current_node: 0,
            freqs,
        }
    }
}
