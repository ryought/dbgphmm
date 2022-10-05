use dbgphmm::common::{ni, Reads};
use dbgphmm::dbg::mocks::mock_intersection_small;
use dbgphmm::vector::{DenseStorage, NodeVec};

fn main() {
    let mut dbg = mock_intersection_small();
    let mut v1: NodeVec<DenseStorage<_>> = NodeVec::new(dbg.n_nodes(), 0.0);
    v1[ni(1)] = 10.0;
    let mut v2: NodeVec<DenseStorage<_>> = NodeVec::new(dbg.n_nodes(), 20.0);
    v2[ni(2)] = 10.0;
    // println!("n={}", dbg.n_nodes());
    // println!("v1={}", v1);
    // println!("v2={}", v2);
    let json = dbg.to_cytoscape_with_attrs_and_historys(
        &[],
        &[],
        &[("v1".to_string(), v1), ("v2".to_string(), v2)],
    );
    println!("{}", json);
}
