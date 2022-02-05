//!
//! Abstraction of types that can be used as an index of vector
//!
pub use petgraph::graph::{EdgeIndex, NodeIndex};

pub trait Indexable: Copy {
    fn new(x: usize) -> Self;
    fn index(&self) -> usize;
}

impl Indexable for usize {
    #[inline]
    fn new(x: usize) -> Self {
        x
    }
    #[inline]
    fn index(&self) -> usize {
        *self
    }
}

impl Indexable for NodeIndex {
    #[inline]
    fn new(x: usize) -> Self {
        NodeIndex::new(x)
    }
    #[inline]
    fn index(&self) -> usize {
        NodeIndex::index(*self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_usize() {
        let x: usize = Indexable::new(10);
        println!("{}", x.index());
    }
    #[test]
    fn index_node_index() {
        let x: NodeIndex = Indexable::new(10);
        println!("{}", x.index());
    }
}
