use super::dbg::{Dbg, DbgEdge, DbgNode};
use super::edge_centric::SimpleEDbg;
use crate::em::extension::flow_intersection::{
    FlowIntersection, FlowIntersectionEdge, FlowIntersectionNode,
};
use crate::graph::Bipartite;
use crate::hmmv2::trans_table::EdgeFreqs;
use crate::kmer::kmer::KmerLike;
use itertools::iproduct;
use petgraph::graph::NodeIndex;

///
/// Intersection information corresponding to a k-1-mer.
///
/// in node-centric de bruijn graph, a k-1-mer can have
/// incoming kmers and outgoing kmers (less than 4 nodes).
///
/// km1mer:    XXX
/// in_nodes:  AXXX, CXXX, GXXX, TXXX
/// out_nodes: XXXA, XXXC, XXXG, XXXT
///
#[derive(Debug, Clone)]
pub struct Intersection<K: KmerLike> {
    bi: Bipartite<K, NodeIndex, ()>,
}

impl<K: KmerLike> Intersection<K> {
    pub fn new(km1mer: K, in_nodes: Vec<NodeIndex>, out_nodes: Vec<NodeIndex>) -> Intersection<K> {
        Intersection {
            bi: Bipartite::from(km1mer, in_nodes, out_nodes, |_, _| ()),
        }
    }
    pub fn km1mer(&self) -> &K {
        &self.bi.id
    }
    pub fn iter_in_nodes(&self) -> impl Iterator<Item = NodeIndex> + '_ {
        self.bi.in_nodes.iter().copied()
    }
    pub fn iter_out_nodes(&self) -> impl Iterator<Item = NodeIndex> + '_ {
        self.bi.out_nodes.iter().copied()
    }
    pub fn in_node(&self, index: usize) -> NodeIndex {
        *self.bi.in_node(index)
    }
    pub fn out_node(&self, index: usize) -> NodeIndex {
        *self.bi.out_node(index)
    }
    pub fn n_in_nodes(&self) -> usize {
        self.bi.n_in()
    }
    pub fn n_out_nodes(&self) -> usize {
        self.bi.n_out()
    }
    pub fn can_uniquely_convertable(&self) -> bool {
        self.n_in_nodes() == 1 || self.n_out_nodes() == 1
    }
    pub fn is_tip_intersection(&self) -> bool {
        self.km1mer().is_null()
    }
}

impl<K: KmerLike> Intersection<K> {
    /// convert a (simple) intersection into the augumented intersection
    /// with flow and freq information.
    pub fn to_flow_intersection<N, E>(
        &self,
        dbg: &Dbg<N, E>,
        freqs: &EdgeFreqs,
    ) -> FlowIntersection<K>
    where
        N: DbgNode,
        E: DbgEdge,
    {
        let in_nodes: Vec<FlowIntersectionNode> = self
            .iter_in_nodes()
            .map(|v| FlowIntersectionNode::new(v, dbg.node(v).copy_num()))
            .collect();
        let out_nodes: Vec<FlowIntersectionNode> = self
            .iter_out_nodes()
            .map(|v| FlowIntersectionNode::new(v, dbg.node(v).copy_num()))
            .collect();
        let edges: Vec<FlowIntersectionEdge> = iproduct!(in_nodes.iter(), out_nodes.iter())
            .map(|(in_node, out_node)| dbg.find_edge(in_node.index, out_node.index).unwrap())
            .map(|e| FlowIntersectionEdge::new(e, freqs[e], dbg.edge(e).copy_num()))
            .collect();
        FlowIntersection::new(self.km1mer().clone(), in_nodes, out_nodes, edges)
    }
}

impl<K: KmerLike> std::fmt::Display for Intersection<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let in_nodes: Vec<NodeIndex> = self.iter_in_nodes().collect();
        let out_nodes: Vec<NodeIndex> = self.iter_out_nodes().collect();
        write!(
            f,
            "Intersection({}) in:{:?} out:{:?}",
            self.km1mer(),
            in_nodes,
            out_nodes
        )
    }
}

///
/// Iterator on de bruijn graph intersections corresponding k-1-mer overlaps
///
/// ## Example (k = 4)
///
/// ATTC ------> TTCG
///
pub struct Intersections<K: KmerLike> {
    edbg: SimpleEDbg<K>,
    current_node: usize,
}

impl<K: KmerLike> Iterator for Intersections<K> {
    type Item = Intersection<K>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.current_node < self.edbg.n_nodes() {
            let intersection = self.edbg.intersection(NodeIndex::new(self.current_node));
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
    pub fn iter_intersections(&self) -> Intersections<N::Kmer> {
        Intersections {
            edbg: self.to_edbg(),
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
            .map(move |i| i.to_flow_intersection(self, freqs))
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
    fn dbg_flow_intersections_simple() {
        let dbg = mock_base();
        let freqs = EdgeFreqs::new(dbg.n_edges(), 1.1);
        println!("{}", dbg);
        println!("{}", freqs);
        for i in dbg.iter_intersections() {
            println!("{}", i);
            let fi = i.to_flow_intersection(&dbg, &freqs);
            println!("{}", fi);

            for (i, j) in iproduct!(0..fi.bi.n_in(), 0..fi.bi.n_out()) {
                let v = fi.bi.in_node(i);
                let w = fi.bi.out_node(j);
                let vw = fi.bi.edge(i, j);
                assert!(dbg.contains_edge(v.index, w.index));
                let e = dbg.find_edge(v.index, w.index).unwrap();
                assert_eq!(freqs[e], vw.freq);
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
            assert!(fi.can_uniquely_convertable());
        }
    }
}
