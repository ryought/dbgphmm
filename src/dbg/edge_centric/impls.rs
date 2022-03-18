use super::{EDbg, EDbgEdge, EDbgNode};
use crate::common::{CopyNum, Sequence};
use crate::kmer::kmer::{Kmer, KmerLike};
use petgraph::graph::DiGraph;
use petgraph::graph::NodeIndex;

/// Basic implementations of EDbg
pub type SimpleEDbg<K> = EDbg<SimpleEDbgNode<K>, SimpleEDbgEdge<K>>;

/// Basic implementations of EDbgNode
pub struct SimpleEDbgNode<K: KmerLike> {
    km1mer: K,
}

impl<K: KmerLike> EDbgNode for SimpleEDbgNode<K> {
    type Kmer = K;
    fn km1mer(&self) -> &K {
        &self.km1mer
    }
}

impl<K: KmerLike> SimpleEDbgNode<K> {
    pub fn new(km1mer: K) -> Self {
        SimpleEDbgNode { km1mer }
    }
}

/// Basic implementations of EDbgNode
pub struct SimpleEDbgEdge<K: KmerLike> {
    kmer: K,
    copy_num: CopyNum,
    origin_node: NodeIndex,
}

impl<K: KmerLike> EDbgEdge for SimpleEDbgEdge<K> {
    type Kmer = K;
    fn kmer(&self) -> &K {
        &self.kmer
    }
    fn copy_num(&self) -> CopyNum {
        self.copy_num
    }
    fn origin_node(&self) -> NodeIndex {
        self.origin_node
    }
}

impl<K: KmerLike> SimpleEDbgEdge<K> {
    pub fn new(kmer: K, copy_num: CopyNum, origin_node: NodeIndex) -> Self {
        SimpleEDbgEdge {
            kmer,
            copy_num,
            origin_node,
        }
    }
}

impl<K: KmerLike> std::fmt::Display for SimpleEDbgNode<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.km1mer)
    }
}

impl<K: KmerLike> std::fmt::Display for SimpleEDbgEdge<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{} (x{}) ({})",
            self.kmer,
            self.copy_num,
            self.origin_node.index()
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
    use crate::kmer::veckmer::VecKmer;

    #[test]
    fn edbg_simple() {
        let mut graph = DiGraph::new();
        let v1 = graph.add_node(SimpleEDbgNode::new(VecKmer::from_bases(b"ATC")));
        let v2 = graph.add_node(SimpleEDbgNode::new(VecKmer::from_bases(b"TCG")));
        let v3 = graph.add_node(SimpleEDbgNode::new(VecKmer::from_bases(b"CGT")));
        let e1 = graph.add_edge(
            v1,
            v2,
            SimpleEDbgEdge::new(VecKmer::from_bases(b"ATCG"), 1, ni(0)),
        );
        let e2 = graph.add_edge(
            v2,
            v3,
            SimpleEDbgEdge::new(VecKmer::from_bases(b"TCGT"), 1, ni(1)),
        );
        let dbg: SimpleEDbg<VecKmer> = SimpleEDbg::new(4, graph);
        println!("{}", dbg);
        assert_eq!(dbg.n_nodes(), 3);
        assert_eq!(dbg.n_edges(), 2);
    }
}
