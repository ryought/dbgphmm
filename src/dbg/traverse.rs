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
        self.traverse_from(NodeIndex::new(0))
    }
    fn traverse_from(&self, from: NodeIndex) -> Traverser<K, N, E> {
        Traverser::new(self, self.to_copy_nums(), from)
    }
    pub fn traverse_all(&self) -> Traveller<K, N, E> {
        Traveller {
            traverser: self.traverse(),
        }
    }
    fn is_valid_cycle(&self, cycle: &Cycle) -> bool {
        let head = cycle.first().unwrap();
        let tail = cycle.last().unwrap();
        self.contains_edge(*tail, *head)
    }
    fn cycle_as_sequence(&self, cycle: &Cycle) -> Sequence {
        assert!(self.is_valid_cycle(cycle));
        cycle.iter().map(|&node| self.emission(node)).collect()
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

/// Iterate through nodes on a cycle from the starting node
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
    fn new(
        dbg: &'a Dbg<K, N, E>,
        copy_nums: NodeCopyNums,
        from: NodeIndex,
    ) -> Traverser<'a, K, N, E> {
        Traverser {
            next_node: Some(from),
            dbg,
            copy_nums,
        }
    }
    fn set_next_node(&mut self, node: NodeIndex) {
        self.next_node = Some(node);
    }
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

pub type Cycle = Vec<NodeIndex>;

pub struct Traveller<'a, K: KmerLike, N: DbgNode<K>, E: DbgEdge> {
    traverser: Traverser<'a, K, N, E>,
}

impl<'a, K, N, E> Traveller<'a, K, N, E>
where
    K: KmerLike,
    N: DbgNode<K>,
    E: DbgEdge,
{
}

impl<'a, K, N, E> Iterator for Traveller<'a, K, N, E>
where
    K: KmerLike,
    N: DbgNode<K>,
    E: DbgEdge,
{
    type Item = Cycle;
    fn next(&mut self) -> Option<Self::Item> {
        // pick a start node
        match self.traverser.find_unvisited_node() {
            Some(node) => {
                // if found, traverse the nodes from the starting node
                self.traverser.set_next_node(node);
                let cycle: Cycle = self.traverser.by_ref().collect();
                // check if this is circular
                Some(cycle)
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
    use crate::common::sequence_to_string;
    use crate::dbg::mocks::*;
    #[test]
    fn dbg_traverse_simple() {
        let dbg = mock_simple();
    }
    #[test]
    fn dbg_traverse_rep() {
        let dbg = mock_rep();
        let circles: Vec<Vec<NodeIndex>> = dbg.traverse_all().collect();
        assert_eq!(circles.len(), 3);
        for circle in circles.iter() {
            println!("{:?}", circle);
            println!("{:?}", sequence_to_string(&dbg.cycle_as_sequence(circle)));
        }
        assert_eq!(dbg.cycle_as_sequence(&circles[0]), b"NNCCCN");
        assert_eq!(dbg.cycle_as_sequence(&circles[1]), b"ANNNAAAAAAAAAAAA");
        assert_eq!(dbg.cycle_as_sequence(&circles[2]), b"CCCCCCCCCCC");
        assert_eq!(circles[0].len() + circles[1].len() + circles[2].len(), 33);
    }
}
