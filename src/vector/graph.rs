//!
//! `Vector` defined on graph
//! Vector that uses NodeIndex or EdgeIndex as an index
//!
//! i.e. `vec[NodeIndex(0)]`
//!
use super::Vector;
pub use petgraph::graph::{EdgeIndex, NodeIndex};

/// Vector that uses petgraph::Node as an index
pub type NodeVec<S> = Vector<S, NodeIndex>;

/// Vector that uses petgraph::Edge as an index
pub type EdgeVec<S> = Vector<S, EdgeIndex>;

//
// converters between min_flow::Flow and EdgeVec<DenseStorage>
//
use crate::min_flow::{Flow, FlowRateLike};
use crate::vector::DenseStorage;

impl<F: FlowRateLike> std::convert::Into<EdgeVec<DenseStorage<F>>> for Flow<F> {
    fn into(self) -> EdgeVec<DenseStorage<F>> {
        EdgeVec::from_inner_vec(self.into())
    }
}

impl<F: FlowRateLike> std::convert::From<EdgeVec<DenseStorage<F>>> for Flow<F> {
    fn from(vector: EdgeVec<DenseStorage<F>>) -> Self {
        vector.to_inner_vec().into()
    }
}

pub fn flow_to_edgevec<F: FlowRateLike>(flow: Flow<F>) -> EdgeVec<DenseStorage<F>> {
    flow.into()
}

pub fn edgevec_to_flow<F: FlowRateLike>(vec: EdgeVec<DenseStorage<F>>) -> Flow<F> {
    vec.into()
}

#[cfg(test)]
mod tests {
    use super::super::dense::DenseStorage;
    use super::*;
    use crate::common::{ei, ni};
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
    #[test]
    fn nodevec_to_edgevec() {
        let mut v: NodeVec<DenseStorage<u32>> = NodeVec::new(5, 0);
        v[ni(1)] = 100;
        let w: EdgeVec<DenseStorage<u32>> = v.switch_index();
        println!("{}", w);
        println!("{}", w[ei(1)]);
        assert_eq!(w[ei(0)], 0);
        assert_eq!(w[ei(1)], 100);
        assert_eq!(w.len(), 5);
    }
}
