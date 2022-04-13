//!
//! traverse related
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use super::intersections::Intersection;
use crate::common::{CopyNum, Sequence};
use crate::kmer::kmer::KmerLike;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

///
/// Path is a list of NodeIndex.
///
pub type Path = Vec<NodeIndex>;

///
/// Path related methods
///
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

///
/// Dbg serializers
///
impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    ///
    /// Convert dbg into sequences
    ///
    pub fn to_seqs(&self) -> Vec<Sequence> {
        self.traverse_all()
            .map(|circle| self.path_as_sequence(&circle))
            .collect()
    }
    ///
    /// Traverse from a single node.
    ///
    pub fn traverse_from(&self, from: NodeIndex) -> Traverser<N, E> {
        Traverser::new(self, self.to_node_copy_nums(), Some(from))
    }
    ///
    /// Traverse all nodes, starting from head nodes.
    ///
    pub fn traverse_all(&self) -> Traveller<N, E> {
        Traveller {
            traverser: Traverser::new(self, self.to_node_copy_nums(), None),
            tips: self.tips(),
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
    fn new(
        dbg: &'a Dbg<N, E>,
        copy_nums: NodeCopyNums,
        next_node: Option<NodeIndex>,
    ) -> Traverser<'a, N, E> {
        Traverser {
            next_node,
            dbg,
            copy_nums,
        }
    }
    fn set_next_node(&mut self, node: NodeIndex) {
        self.next_node = Some(node);
    }
    fn is_unvisited(&self, node: NodeIndex) -> bool {
        self.copy_nums[node] > 0
    }
    fn find_unvisited_child(&self, node: NodeIndex) -> Option<NodeIndex> {
        self.dbg
            .childs(node)
            .map(|(_, child, _)| child)
            .find(|&child| self.is_unvisited(child))
    }
    fn find_unvisited_node(&self) -> Option<NodeIndex> {
        self.dbg
            .nodes()
            .map(|(node, _)| node)
            .find(|&node| self.is_unvisited(node))
    }
    fn as_path(&mut self) -> Path {
        self.collect()
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

/// Eulerian Traverse all nodes in dbg from starting nodes.
///
/// * as iterator
/// * as vector of paths by using `as_paths()`
///
pub struct Traveller<'a, N: DbgNode, E: DbgEdge> {
    traverser: Traverser<'a, N, E>,
    tips: Intersection<N::Kmer>,
}

impl<'a, N, E> Traveller<'a, N, E>
where
    N: DbgNode,
    E: DbgEdge,
{
    /// Find the starting node
    ///
    /// * Use head nodes if unvisited.
    /// * Otherwise, use unvisited ordinary node
    fn find_starting_node(&self) -> Option<NodeIndex> {
        // (1) find unvisited head node (e.g. NNNA)
        for head_node in self.tips.iter_out_nodes() {
            if self.traverser.is_unvisited(head_node) {
                return Some(head_node);
            }
        }

        // (2) if head node is not found, then use unvisited ordinary node
        self.traverser.find_unvisited_node()
    }
    /// Get paths
    ///
    pub fn as_paths(&mut self) -> Vec<Path> {
        self.collect()
    }
}

impl<'a, N, E> Iterator for Traveller<'a, N, E>
where
    N: DbgNode,
    E: DbgEdge,
{
    type Item = Path;
    fn next(&mut self) -> Option<Self::Item> {
        // pick a start node
        match self.find_starting_node() {
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
    use crate::common::{ni, sequence_to_string, Reads};
    use crate::dbg::impls::SimpleDbg;
    use crate::dbg::mocks::*;
    use crate::kmer::VecKmer;
    #[test]
    fn dbg_traverse_simple() {
        let dbg = mock_simple();
        println!("{}", dbg);

        let p = dbg.traverse_from(ni(0)).as_path();
        println!("{:?}", p);
        assert_eq!(
            p,
            vec![
                ni(0),
                ni(1),
                ni(10),
                ni(9),
                ni(8),
                ni(4),
                ni(11),
                ni(2),
                ni(3),
                ni(5),
                ni(7),
                ni(12),
                ni(13),
                ni(6)
            ]
        );

        let v = dbg.traverse_all().find_starting_node();
        assert!(v.is_some());
        assert_eq!(v.unwrap(), ni(10));
        assert_eq!(dbg.kmer(v.unwrap()), &VecKmer::from_bases(b"NNNA"));

        let ps = dbg.traverse_all().as_paths();
        println!("{:?}", ps);
        assert_eq!(
            ps,
            vec![vec![
                ni(10),
                ni(9),
                ni(8),
                ni(4),
                ni(11),
                ni(2),
                ni(3),
                ni(5),
                ni(7),
                ni(12),
                ni(13),
                ni(6),
                ni(0),
                ni(1),
            ]]
        );

        let seqs = dbg.to_seqs();
        for seq in seqs.iter() {
            println!("{}", sequence_to_string(seq));
        }
        assert_eq!(seqs, vec![b"AAAGCTTGATTNNN"]);
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
        assert_eq!(dbg.path_as_sequence(&circles[0]), b"CCCNNNAAANNN");
        assert_eq!(dbg.path_as_sequence(&circles[1]), b"AAAAAAAAAA");
        assert_eq!(dbg.path_as_sequence(&circles[2]), b"CCCCCCCCCCC");
        assert_eq!(circles[0].len() + circles[1].len() + circles[2].len(), 33);
    }
    #[test]
    fn dbg_traverse_to_seqs() {
        let dbg: SimpleDbg<VecKmer> =
            SimpleDbg::from_reads(8, &Reads::from(vec![b"ATTCGATCGAT".to_vec()]));
        println!("{}", dbg);
        for seq in dbg.to_seqs().iter() {
            println!("{}", sequence_to_string(seq));
        }
    }
}
