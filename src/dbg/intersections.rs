use super::dbg::{Dbg, DbgEdge, DbgNode};
use super::edge_centric::{SimpleEDbg, SimpleEDbgEdge, SimpleEDbgNode};
use crate::kmer::kmer::KmerLike;
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
    pub km1mer: K,
    pub in_nodes: Vec<NodeIndex>,
    pub out_nodes: Vec<NodeIndex>,
}

impl<K: KmerLike> Intersection<K> {
    pub fn n_in_nodes(&self) -> usize {
        self.in_nodes.len()
    }
    pub fn n_out_nodes(&self) -> usize {
        self.out_nodes.len()
    }
}

impl<K: KmerLike> std::fmt::Display for Intersection<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "Intersection({}) in:{:?} out:{:?}",
            self.km1mer, self.in_nodes, self.out_nodes
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
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ni;
    use crate::dbg::hashdbg_v2::HashDbg;
    use crate::dbg::impls::{SimpleDbg, SimpleDbgEdge, SimpleDbgNode};
    use crate::kmer::veckmer::VecKmer;

    #[test]
    fn dbg_intersections() {
        let hd: HashDbg<VecKmer> = HashDbg::from_seq(4, b"ATCGGCT");
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_hashdbg(&hd);
        println!("{}", dbg);
        for i in dbg.iter_intersections() {
            for &in_node in i.in_nodes.iter() {
                assert_eq!(dbg.node(in_node).kmer().suffix(), i.km1mer);
            }
            for &out_node in i.out_nodes.iter() {
                assert_eq!(dbg.node(out_node).kmer().prefix(), i.km1mer);
            }
            assert_eq!(i.n_in_nodes(), 1);
            assert_eq!(i.n_out_nodes(), 1);
            println!("{}", i);
        }
    }
}
