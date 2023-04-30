use crate::utils::{check_memory_usage, timer};
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};

trait A: std::marker::Sized {
    fn a<N: A, E: A>(&self, graph: &DiGraph<N, E>) -> usize;
}

trait GA {
    fn ga(&self) -> usize;
}

impl<N: A, E: A> GA for DiGraph<N, E> {
    fn ga(&self) -> usize {
        let n = self.node_weight(NodeIndex::new(0)).unwrap().a(self);
        let e = self.edge_weight(EdgeIndex::new(0)).unwrap().a(self);
        n + e
    }
}

fn t<N: A, E: A>(g: DiGraph<N, E>) {}

// (1)
trait PN {
    fn init_prob<N: PN, E: PE>(&self, node: NodeIndex, graph: &DiGraph<N, E>) -> f64;
}
trait PE {
    fn trans_prob<N: PN, E: PE>(&self, edge: EdgeIndex, graph: &DiGraph<N, E>) -> f64;
}
fn to_frozen_type1<N, E>(graph: &DiGraph<N, E>)
where
    N: PN,
    E: PE,
{
}
// (2)
trait PG {
    type N;
    type E;
    fn graph(&self) -> &DiGraph<Self::N, Self::E>;
    fn init_prob(&self, node: NodeIndex) -> f64;
    fn trans_prob(&self, edge: EdgeIndex) -> f64;
}
impl<N: A, E: A> PG for DiGraph<N, E> {
    type N = N;
    type E = E;
    fn graph(&self) -> &DiGraph<N, E> {
        self
    }
    fn init_prob(&self, node: NodeIndex) -> f64 {
        1.0
    }
    fn trans_prob(&self, edge: EdgeIndex) -> f64 {
        0.0
    }
}
fn to_frozen_type2<G: PG>(model: G) {
    // let g2 = model.graph().map();
}

/*
// you can use Self below by restricting trait instances as a member of Copy or Sized
pub trait PHMMLikeNode: std::marker::Sized {
    fn base(&self) -> u8;
    // fn init_prob<E>(&self, graph: &DiGraph<Self, E>) -> Prob;
    fn init_prob(&self) -> Prob;
}
pub trait PHMMLikeEdge: std::marker::Sized {
    // fn trans_prob<N>(&self, graph: &DiGraph<N, Self>) -> Prob;
    fn trans_prob(&self) -> Prob;
}
pub trait PHMMLikeGraph {
    fn to_phmm(&self, param: PHMMParams) -> PModel;
}
impl<N: PHMMLikeNode, E: PHMMLikeEdge> PHMMLikeGraph for DiGraph<N, E> {
    fn to_phmm(&self, param: PHMMParams) -> PModel {
        let graph = self.map(
            |_, vw| PNode::new(0, vw.init_prob(), true, vw.base()),
            |_, ew| PEdge::new(ew.trans_prob()),
        );
        PModel { param, graph }
    }
}
// for Dbg
impl<N: DbgNode, E: DbgEdge> PHMMLikeGraph for Dbg<N, E> {
    fn to_phmm(&self, param: PHMMParams) -> PModel {
        self.graph.to_phmm(param)
    }
}
impl<N: DbgNode> PHMMLikeNode for N {
    fn base(&self) -> u8 {
        self.emission()
    }
    fn init_prob(&self) -> Prob {
        // TODO
        Prob::from_prob(0.0)
    }
}
impl<E: DbgEdge> PHMMLikeEdge for E {
    fn trans_prob(&self) -> Prob {
        // TODO
        Prob::from_prob(0.0)
    }
}
*/

fn p<I>(i: I)
where
    I: IntoIterator,
    I::Item: AsRef<(usize, u8)>,
{
    for item in i {
        let (x, y) = item.as_ref();
        println!("x={} y={}", x, y);
    }
}
fn f() {
    let v: Vec<(usize, u8)> = vec![(10, 2), (11, 5)];
    // this does not compile
    // p(&v);
}

fn create_large_graph() -> DiGraph<usize, ()> {
    // create graph with 10M nodes
    let mut g = DiGraph::new();
    for i in 0..1_000_000_000 {
        g.add_node(i);
    }
    g
}

fn test_memory_usage_of_large_graph() {
    let (g, t) = timer(|| {
        check_memory_usage();
        let g = create_large_graph();
        check_memory_usage();
        g
    });
    println!("{}", t);
    println!("{}", g.node_count());
}

//
// test matrix determinant using openblas and ndarray
//
use ndarray::*;
use ndarray_linalg::*;
fn calc_det() {
    let a = arr2(&[[2.0, 1.0, 2.0], [-2.0, 2.0, 1.0], [1.0, 2.0, -2.0]]);
    println!("a={}", a);
    let d = a.det().unwrap();
    println!("det(a)={}", d);
    assert_eq!(d, -27.);
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn calc_det_test() {
        calc_det();
    }
}
