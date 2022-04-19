use super::dbg::{Dbg, DbgEdge, DbgNode};
use super::edge_centric::{IntersectionBase, SimpleEDbg};
use crate::dbg::flow_intersection::{FlowIntersection, FlowIntersectionEdge, FlowIntersectionNode};
use crate::graph::Bipartite;
use crate::hmmv2::trans_table::EdgeFreqs;
use crate::kmer::kmer::KmerLike;
use itertools::iproduct;
use petgraph::graph::NodeIndex;

///
/// Iterator on de bruijn graph intersections corresponding k-1-mer overlaps
///
/// ## Example (k = 4)
///
/// ATTC ------> TTCG
///
pub struct Intersections<'a, N: DbgNode, E: DbgEdge> {
    edbg: SimpleEDbg<N::Kmer>,
    dbg: &'a Dbg<N, E>,
    current_node: usize,
}

impl<'a, N: DbgNode, E: DbgEdge> Iterator for Intersections<'a, N, E> {
    type Item = FlowIntersection<N::Kmer>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.current_node < self.edbg.n_nodes() {
            let intersection_base = self.edbg.intersection(NodeIndex::new(self.current_node));
            let intersection = self.dbg.to_flow_intersection(&intersection_base, None);
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
    pub fn iter_intersections(&self) -> Intersections<N, E> {
        Intersections {
            edbg: self.to_edbg(),
            dbg: self,
            current_node: 0,
        }
    }
    ///
    /// get an iterator over intersections (k-1-mer)
    /// **with flow and copy_nums information**
    ///
    pub fn iter_flow_intersections<'a>(
        &'a self,
        freqs: &'a EdgeFreqs,
    ) -> impl Iterator<Item = FlowIntersection<N::Kmer>> + 'a {
        self.iter_intersections()
            .map(move |i| i.augment_freqs(freqs))
    }
    /// convert a (simple) intersection into the augumented intersection
    /// with flow and freq information.
    pub fn to_flow_intersection(
        &self,
        intersection_base: &IntersectionBase<N::Kmer>,
        freqs: Option<&EdgeFreqs>,
    ) -> FlowIntersection<N::Kmer> {
        let in_nodes: Vec<FlowIntersectionNode> = intersection_base
            .in_nodes()
            .iter()
            .map(|&v| FlowIntersectionNode::new(v, self.node(v).copy_num()))
            .collect();
        let out_nodes: Vec<FlowIntersectionNode> = intersection_base
            .out_nodes()
            .iter()
            .map(|&v| FlowIntersectionNode::new(v, self.node(v).copy_num()))
            .collect();
        let edges: Vec<FlowIntersectionEdge> = iproduct!(in_nodes.iter(), out_nodes.iter())
            .map(|(in_node, out_node)| self.find_edge(in_node.index, out_node.index).unwrap())
            .map(|e| {
                let freq = match freqs {
                    Some(freqs) => Some(freqs[e]),
                    None => None,
                };
                FlowIntersectionEdge::new(e, freq, self.edge(e).copy_num())
            })
            .collect();
        FlowIntersection::new(
            intersection_base.km1mer().clone(),
            in_nodes,
            out_nodes,
            edges,
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
    use crate::dbg::mocks::{mock_base, mock_intersection};
    use crate::kmer::veckmer::VecKmer;

    #[test]
    fn dbg_intersections_simple() {
        let dbg = mock_base();
        println!("{}", dbg);
        for i in dbg.iter_intersections() {
            for in_node in i.iter_in_nodes() {
                assert_eq!(&dbg.node(in_node).kmer().suffix(), i.km1mer());
            }
            for out_node in i.iter_out_nodes() {
                assert_eq!(&dbg.node(out_node).kmer().prefix(), i.km1mer());
            }
            assert_eq!(i.n_in_nodes(), 1);
            assert_eq!(i.n_out_nodes(), 1);
            println!("{}", i);
        }
    }
    #[test]
    fn dbg_intersections_twoseqs() {
        let dbg = mock_intersection();
        for i in dbg.iter_intersections() {
            if i.km1mer() == &VecKmer::from_bases(b"nnn") {
                assert_eq!(i.n_in_nodes(), 2);
                assert_eq!(i.n_out_nodes(), 2);
            } else if i.km1mer() == &VecKmer::from_bases(b"TAG") {
                assert_eq!(i.n_in_nodes(), 2);
                assert_eq!(i.n_out_nodes(), 2);
            } else {
                assert_eq!(i.n_in_nodes(), 1);
                assert_eq!(i.n_out_nodes(), 1);
            }

            for in_node in i.iter_in_nodes() {
                assert_eq!(&dbg.node(in_node).kmer().suffix(), i.km1mer());
            }
            for out_node in i.iter_out_nodes() {
                assert_eq!(&dbg.node(out_node).kmer().prefix(), i.km1mer());
            }
            println!("{}", i);
        }
    }
    #[test]
    fn dbg_flow_intersections_simple_0() {
        let dbg = mock_base();
        let freqs = EdgeFreqs::new(dbg.n_edges(), 1.1);
        println!("{}", dbg);
        println!("{}", freqs);
        for i in dbg.iter_intersections() {
            println!("without-flow {}", i);
            let fi = i.augment_freqs(&freqs);
            println!("with-flow {}", fi);

            for (i, j) in iproduct!(0..fi.bi.n_in(), 0..fi.bi.n_out()) {
                let v = fi.bi.in_node(i);
                let w = fi.bi.out_node(j);
                let vw = fi.bi.edge(i, j);
                assert!(dbg.contains_edge(v.index, w.index));
                let e = dbg.find_edge(v.index, w.index).unwrap();
                assert_eq!(Some(freqs[e]), vw.freq);
            }
        }
    }
    #[test]
    fn dbg_flow_intersections_simple_iter() {
        let dbg = mock_base();
        let freqs = EdgeFreqs::new(dbg.n_edges(), 1.1);
        println!("{}", dbg);
        println!("{}", freqs);
        for fi in dbg.iter_flow_intersections(&freqs) {
            println!("{}", fi);
            assert_eq!(fi.n_in_nodes(), 1);
            assert_eq!(fi.n_out_nodes(), 1);
            assert!(fi.has_valid_node_copy_nums());
            assert!(fi.can_unique_resolvable());
        }
    }
}
