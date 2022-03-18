//!
//! traverse related
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use crate::common::{CopyNum, Sequence};
use crate::kmer::kmer::KmerLike;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

///
/// list of NodeIndex
///
pub type Path = Vec<NodeIndex>;

/// Path related methods
impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    ///
    /// check if the given path (node list) is cycle.
    ///
    pub fn is_cycle(&self, path: &Path) -> bool {
        let head = path.first().unwrap();
        let tail = path.last().unwrap();
        self.contains_edge(*tail, *head)
    }
    ///
    /// convert node list into bases
    ///
    pub fn path_as_sequence(&self, path: &Path) -> Sequence {
        path.iter().map(|&node| self.emission(node)).collect()
    }
}

impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    /// Convert dbg into sequences
    ///
    pub fn to_seqs(&self) -> Vec<Sequence> {
        self.traverse_all()
            .map(|circle| self.path_as_sequence(&circle))
            .collect()
    }
    /// Traverse all nodes
    ///
    /// ## TODO
    ///
    /// * always should start from head nodes.
    ///
    pub fn traverse(&self) -> Traverser<N, E> {
        self.traverse_from(NodeIndex::new(0))
    }
    fn traverse_from(&self, from: NodeIndex) -> Traverser<N, E> {
        Traverser::new(self, self.to_node_copy_nums(), from)
    }
    pub fn traverse_all(&self) -> Traveller<N, E> {
        Traveller {
            traverser: self.traverse(),
        }
    }
}

/// Iterate through nodes on a cycle from the starting node
pub struct Traverser<'a, N: DbgNode, E: DbgEdge> {
    /// the node index it will visit in the next step
    next_node: Option<NodeIndex>,
    /// reference to the dbg struct
    dbg: &'a Dbg<N, E>,
    /// remaining unvisited copy numbers of each nodes
    copy_nums: NodeCopyNums,
}

impl<'a, N, E> Traverser<'a, N, E>
where
    N: DbgNode,
    E: DbgEdge,
{
    fn new(dbg: &'a Dbg<N, E>, copy_nums: NodeCopyNums, from: NodeIndex) -> Traverser<'a, N, E> {
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

impl<'a, N, E> Iterator for Traverser<'a, N, E>
where
    N: DbgNode,
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

pub struct Traveller<'a, N: DbgNode, E: DbgEdge> {
    traverser: Traverser<'a, N, E>,
}

impl<'a, N, E> Traveller<'a, N, E>
where
    N: DbgNode,
    E: DbgEdge,
{
}

impl<'a, N, E> Iterator for Traveller<'a, N, E>
where
    N: DbgNode,
    E: DbgEdge,
{
    type Item = Path;
    fn next(&mut self) -> Option<Self::Item> {
        // pick a start node
        match self.traverser.find_unvisited_node() {
            Some(node) => {
                // if found, traverse the nodes from the starting node
                self.traverser.set_next_node(node);
                let cycle: Path = self.traverser.by_ref().collect();
                // check if this is circular
                assert!(self.traverser.dbg.is_cycle(&cycle));
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
        let seqs = dbg.to_seqs();
        println!("{}", dbg);
        for seq in seqs.iter() {
            println!("{}", sequence_to_string(seq));
        }
        assert_eq!(seqs, vec![b"NNAAAGCTTGATTN"]);
    }
    #[test]
    fn dbg_traverse_rep() {
        let dbg = mock_rep();
        let circles: Vec<Vec<NodeIndex>> = dbg.traverse_all().collect();
        assert_eq!(circles.len(), 3);
        for circle in circles.iter() {
            println!("{:?}", circle);
            println!("{:?}", sequence_to_string(&dbg.path_as_sequence(circle)));
        }
        assert_eq!(dbg.path_as_sequence(&circles[0]), b"NNCCCN");
        assert_eq!(dbg.path_as_sequence(&circles[1]), b"ANNNAAAAAAAAAAAA");
        assert_eq!(dbg.path_as_sequence(&circles[2]), b"CCCCCCCCCCC");
        assert_eq!(circles[0].len() + circles[1].len() + circles[2].len(), 33);
    }
}
