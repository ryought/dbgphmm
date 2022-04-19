//!
//! traverse related
//!
use super::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use crate::common::{CopyNum, SeqStyle, Sequence, StyledSequence, NULL_BASE};
use crate::dbg::flow_intersection::FlowIntersection;
use crate::kmer::kmer::KmerLike;
use itertools::Itertools;
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
    /// the given path is linear or not, that is
    /// starts from NNNX and ends with YNNN.
    ///
    pub fn is_linear(&self, path: &Path) -> bool {
        let head = path.first().unwrap();
        let tail = path.last().unwrap();
        self.kmer(*head).is_head() && self.kmer(*tail).is_tail()
    }
    ///
    /// convert node list into bases
    ///
    pub fn path_as_sequence(&self, path: &Path) -> Sequence {
        path.iter().map(|&node| self.emission(node)).collect()
    }
    ///
    /// convert node list into bases
    /// with omitting null emissions
    ///
    pub fn path_as_sequence_without_null(&self, path: &Path) -> Sequence {
        path.iter()
            .map(|&node| self.emission(node))
            .filter(|&base| base != NULL_BASE)
            .collect()
    }
    ///
    /// convert a path (node list) into sequence (vec of bases) with style.
    ///
    /// the style can be either linear or circular
    ///
    pub fn path_as_styled_sequence(&self, path: &Path) -> StyledSequence {
        let style = if self.is_linear(path) {
            SeqStyle::Linear
        } else {
            SeqStyle::Circular
        };
        StyledSequence::new(self.path_as_sequence_without_null(path), style)
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
    /// Convert dbg into styled sequences
    ///
    pub fn to_styled_seqs(&self) -> Vec<StyledSequence> {
        let (styled_seqs, _) = self.to_styled_seqs_with_n_choices();
        styled_seqs
    }
    ///
    /// Convert dbg into styled sequences, with n_choices information.
    ///
    pub fn to_styled_seqs_with_n_choices(&self) -> (Vec<StyledSequence>, usize) {
        let mut traveller = self.traverse_all();
        let styled_seqs: Vec<StyledSequence> = traveller
            .by_ref()
            .map(|circle| self.path_as_styled_sequence(&circle))
            .collect();
        let n_choices = traveller.n_path_choices();
        (styled_seqs, n_choices)
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
    /// the number of choices in the next step
    n_next_node_choices: usize,
    /// reference to the dbg struct
    dbg: &'a Dbg<N, E>,
    /// remaining unvisited copy numbers of each nodes
    copy_nums: NodeCopyNums,
    /// the number of alternative paths to current node.
    n_path_choices: usize,
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
            n_next_node_choices: 1,
            dbg,
            copy_nums,
            n_path_choices: 1,
        }
    }
    fn set_next_node(&mut self, node: NodeIndex) {
        self.next_node = Some(node);
        self.n_next_node_choices = 1;
    }
    fn is_unvisited(&self, node: NodeIndex) -> bool {
        self.copy_nums[node] > 0
    }
    ///
    /// Find the unvisited (= copy_nums of the node is remaining) child
    /// of the node and the number of choices of them.
    ///
    /// to achieve stable traversing (i.e. the same traversal order even though the
    /// index of nodes/edges are different), the child indexes is returned in order of
    /// k-mer dictionary order.
    ///
    fn find_unvisited_child(&self, node: NodeIndex) -> (Option<NodeIndex>, usize) {
        if self.dbg.kmer(node).is_tail() {
            // if the node is tail (XNNNN), break the path.
            (None, 0)
        } else {
            let mut unvisited_childs: Vec<_> = self
                .dbg
                .childs(node)
                .map(|(_, child, _)| child)
                .filter(|&child| self.is_unvisited(child))
                .collect();
            unvisited_childs.sort_by_key(|&v| self.dbg.kmer(v));
            (unvisited_childs.first().copied(), unvisited_childs.len())
        }
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
                self.n_path_choices = self.n_path_choices.saturating_mul(self.n_next_node_choices);

                // search for new next node
                // and store the n_choices.
                let (next_node, n_choices) = self.find_unvisited_child(node);
                self.next_node = next_node;
                self.n_next_node_choices = n_choices;

                Some(node)
            }
            None => None,
        }
    }
}

/// Eulerian Traverse all nodes in dbg from starting nodes.
///
/// It can (only) be created from dbg using `Dbg::traverse_all()`.
///
/// * as iterator
/// * as vector of paths by using `as_paths()`
///
pub struct Traveller<'a, N: DbgNode, E: DbgEdge> {
    traverser: Traverser<'a, N, E>,
    tips: FlowIntersection<N::Kmer>,
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
    /// Access to the current n_path_choices
    /// in the traverser.
    pub fn n_path_choices(&self) -> usize {
        self.traverser.n_path_choices
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

        // #1 traverse starting from ni(0): ni(1) is tail
        let p = dbg.traverse_from(ni(0)).as_path();
        println!("{:?}", p);
        println!("{}", dbg.kmer(ni(0)));
        println!("{}", dbg.kmer(ni(1)));
        assert_eq!(p, vec![ni(0), ni(1)]);

        // #2 traverse starting from head
        let v = dbg.traverse_all().find_starting_node();
        assert!(v.is_some());
        assert_eq!(v.unwrap(), ni(10));
        assert_eq!(dbg.kmer(v.unwrap()), &VecKmer::from_bases(b"nnnA"));

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
        assert_eq!(seqs, vec![b"AAAGCTTGATTnnn"]);

        let seqs = dbg.to_styled_seqs();
        for seq in seqs.iter() {
            println!("{}", seq);
        }
        assert_eq!(seqs.len(), 1);
        assert_eq!(format!("{}", seqs[0]), "L:AAAGCTTGATT");

        println!("{}", dbg);
        assert_eq!(dbg.to_string(), "4,L:AAAGCTTGATT");
    }
    #[test]
    fn dbg_traverse_rep() {
        let dbg = mock_rep();
        println!("{}", dbg.to_dot());

        let circles: Vec<Vec<NodeIndex>> = dbg.traverse_all().collect();
        assert_eq!(circles.len(), 2);
        for circle in circles.iter() {
            println!("{:?}", circle);
            println!("{:?}", sequence_to_string(&dbg.path_as_sequence(circle)));
        }
        assert_eq!(dbg.path_as_sequence(&circles[0]), b"CCCCCCCCCCCCCCnnn");
        assert_eq!(dbg.path_as_sequence(&circles[1]), b"AAAAAAAAAAAAAnnn");
        assert_eq!(circles[0].len() + circles[1].len(), 33);

        let seqs = dbg.to_styled_seqs();
        for seq in seqs.iter() {
            println!("{}", seq);
        }
        assert_eq!(seqs.len(), 2);
        assert_eq!(format!("{}", seqs[0]), "L:CCCCCCCCCCCCCC");
        assert_eq!(format!("{}", seqs[1]), "L:AAAAAAAAAAAAA");

        println!("{}", dbg);
        assert_eq!(dbg.to_string(), "4,L:CCCCCCCCCCCCCC,L:AAAAAAAAAAAAA");
    }
    #[test]
    fn dbg_traverse_to_seqs() {
        let dbg: SimpleDbg<VecKmer> =
            SimpleDbg::from_reads(8, &Reads::from(vec![b"ATTCGATCGAT".to_vec()]));
        println!("{}", dbg);
        for seq in dbg.to_seqs().iter() {
            println!("{}", sequence_to_string(seq));
        }

        let seqs = dbg.to_styled_seqs();
        for seq in seqs.iter() {
            println!("{}", seq);
        }
        assert_eq!(seqs.len(), 1);
        assert_eq!(format!("{}", seqs[0]), "L:ATTCGATCGAT");
    }
    #[test]
    fn dbg_traverse_n_choices() {
        for (dbg, n_choices_true) in [
            (mock_simple(), 1),
            (mock_rep(), 2097152),
            (mock_intersection(), 2),
        ] {
            let (seqs, n_choices) = dbg.to_styled_seqs_with_n_choices();
            for seq in seqs.iter() {
                println!("{}", seq);
            }
            println!("n_choices={}", n_choices);
            assert_eq!(n_choices, n_choices_true);
        }
    }
}
