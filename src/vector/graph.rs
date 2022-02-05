//!
//! `Vector` defined on graph
//! Vector that uses NodeIndex or EdgeIndex as an index
//!
//! i.e. `vec[NodeIndex(0)]`
//!
use super::{Storage, Vector};
pub use petgraph::graph::{EdgeIndex, NodeIndex};
use std::ops::{Add, AddAssign, Index, IndexMut};

///
/// Vector that supports index access by petgraph::NodeIndex
///
#[derive(Debug, Clone)]
pub struct NodeVec<S: Storage>(pub Vector<S>);

impl<S: Storage> Index<NodeIndex> for NodeVec<S> {
    type Output = S::Item;
    fn index(&self, node: NodeIndex) -> &Self::Output {
        &self.0[node.index()]
    }
}

impl<S: Storage> IndexMut<NodeIndex> for NodeVec<S> {
    fn index_mut(&mut self, node: NodeIndex) -> &mut Self::Output {
        &mut self.0[node.index()]
    }
}

impl<S: Storage> NodeVec<S> {
    /// Create a new Vector, with fixed size and filled by default_value.
    pub fn new(size: usize, default_value: S::Item) -> NodeVec<S> {
        NodeVec(Vector::new(size, default_value))
    }
    /// get the size of the vector
    pub fn len(&self) -> usize {
        self.0.len()
    }
}

impl<S: Storage> NodeVec<S> {
    /// Get an iterator on (index, item).
    pub fn iter<'a>(&'a self) -> impl 'a + Iterator<Item = (NodeIndex, S::Item)> {
        self.0.iter().map(|(i, v)| (NodeIndex::new(i), v))
    }
}

/*
impl<'a, 'b, S> Add<&'a NodeVec<S>> for &'b NodeVec<S>
where
    S: Storage,
    S::Item: Add<Output = S::Item>,
{
    type Output = NodeVec<S>;
    fn add(self, other: &'a NodeVec<S>) -> Self::Output {
        let v = &self.0 + &other.0;
        NodeVec(v)
    }
}

impl<'a, S> AddAssign<&'a NodeVec<S>> for NodeVec<S>
where
    S: IterableStorage<'a>,
    S::Item: Add<Output = S::Item>,
{
    fn add_assign(&mut self, other: &'a NodeVec<S>) {
        self.0 += &other.0;
    }
}
*/

#[cfg(test)]
mod tests {
    use super::super::dense::DenseStorage;
    use super::*;
    use crate::prob::Prob;

    #[test]
    fn nodevec() {
        let mut v: NodeVec<DenseStorage<u32>> = NodeVec::new(5, 0);
        v[NodeIndex::new(1)] = 100;
        assert_eq!(v.len(), 5);
        assert_eq!(v[NodeIndex::new(0)], 0);
        assert_eq!(v[NodeIndex::new(1)], 100);
        let w: Vec<(NodeIndex, u32)> = v.iter().collect();
        assert_eq!(
            w,
            vec![
                (NodeIndex::new(0), 0),
                (NodeIndex::new(1), 100),
                (NodeIndex::new(2), 0),
                (NodeIndex::new(3), 0),
                (NodeIndex::new(4), 0),
            ]
        );

        let mut v2: NodeVec<DenseStorage<u32>> = NodeVec::new(5, 0);
        v2[NodeIndex::new(0)] = 222;

        let added = &v + &v2;
        assert_eq!(added[NodeIndex::new(0)], 0 + 222);
        assert_eq!(added[NodeIndex::new(1)], 100);
        assert_eq!(added[NodeIndex::new(2)], 0);
        assert_eq!(added[NodeIndex::new(3)], 0);
    }
    #[test]
    fn nodevec_prob() {
        let mut v1: NodeVec<DenseStorage<Prob>> = NodeVec::new(5, Prob::from_prob(0.0));
        v1[NodeIndex::new(1)] = Prob::from_prob(0.1);
        v1[NodeIndex::new(3)] = Prob::from_prob(0.2);
        println!("{:?}", v1);
        let mut v2: NodeVec<DenseStorage<Prob>> = NodeVec::new(5, Prob::from_prob(0.0));
        v2[NodeIndex::new(1)] = Prob::from_prob(0.2);
        v2[NodeIndex::new(4)] = Prob::from_prob(0.3);
        // v1 += &v2;
        println!("{:?}", v1);
    }
}
