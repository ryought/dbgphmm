//!
//! traverse related
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use crate::common::{CopyNum, Sequence};
use crate::kmer::kmer::KmerLike;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

impl<K: KmerLike, N: DbgNode<K>, E: DbgEdge> Dbg<K, N, E> {
    /// TODO
    pub fn to_seqs(&self) -> Vec<Sequence> {
        unimplemented!();
    }
    ///
    ///
    ///
    pub fn traverse(&self) -> Traverser<K, N, E> {
        Traverser {
            // should start from root
            next_node: Some(NodeIndex::new(0)),
            dbg: self,
            copy_nums: self.to_copy_nums(),
        }
    }
    ///
    /// Create the NodeVec with copy numbers of each node
    ///
    pub fn to_copy_nums(&self) -> NodeCopyNums {
        let mut v: NodeCopyNums = NodeCopyNums::new(self.n_nodes(), 0);
        for (node, weight) in self.nodes() {
            v[node] = weight.copy_num();
        }
        v
    }
}

pub struct Traverser<'a, K: KmerLike, N: DbgNode<K>, E: DbgEdge> {
    /// the node index it will visit in the next step
    next_node: Option<NodeIndex>,
    /// reference to the dbg struct
    dbg: &'a Dbg<K, N, E>,
    /// remaining unvisited copy numbers of each nodes
    copy_nums: NodeCopyNums,
}

impl<'a, K, N, E> Traverser<'a, K, N, E>
where
    K: KmerLike,
    N: DbgNode<K>,
    E: DbgEdge,
{
    fn find_unvisited_child(&self, node: NodeIndex) -> Option<NodeIndex> {
        self.dbg
            .childs(node)
            .map(|(_, child, _)| child)
            .find(|child| self.copy_nums[*child] > 0)
    }
    fn find_unvisited_node(&self) -> Option<NodeIndex> {
        self.dbg
            .nodes()
            .map(|(node, _)| node)
            .find(|node| self.copy_nums[*node] > 0)
    }
}

impl<'a, K, N, E> Iterator for Traverser<'a, K, N, E>
where
    K: KmerLike,
    N: DbgNode<K>,
    E: DbgEdge,
{
    type Item = NodeIndex;
    fn next(&mut self) -> Option<Self::Item> {
        match self.next_node {
            Some(node) => {
                // visit this node
                self.copy_nums[node] -= 1;
                // search for new next node
                self.next_node = self.find_unvisited_child(node);
                Some(node)
            }
            None => None,
        }
    }
}

//
// tests
//
#[cfg(test)]
mod tests {
    use super::*;
    use crate::dbg::hashdbg_v2::HashDbg;
    use crate::dbg::impls::SimpleDbg;
    use crate::kmer::veckmer::VecKmer;
    #[test]
    fn dbg_traverse() {
        let hd: HashDbg<VecKmer> = HashDbg::from_seq(4, b"AAAGCTTGATT");
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_hashdbg(&hd);
        println!("{}", dbg);
        for node in dbg.traverse() {
            println!("{:?}", node);
        }
    }
}
