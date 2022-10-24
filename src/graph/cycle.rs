//!
//! Cycle in graph
//!
use crate::dbg::dbg::{EdgeCopyNums, NodeCopyNums};
use fixedbitset::FixedBitSet;
use itertools::Itertools;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex, UnGraph};
use petgraph::visit::EdgeRef;
use std::cmp::Ordering;

///
/// cycle (as a list of edges)
///
#[derive(Debug, Clone, PartialEq)]
pub struct Cycle(Vec<EdgeIndex>);
//
// Cycle
//
impl std::fmt::Display for Cycle {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.0.iter().map(|e| e.index()).join(","))
    }
}

///
/// simple cycle (without edge repetition)
///
#[derive(Debug, Clone, PartialEq)]
pub struct SimpleCycle(FixedBitSet);
//
// TODO
// should use to_cycle and show as edge list.
//
impl std::fmt::Display for SimpleCycle {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[")?;
        for i in 0..self.0.len() {
            if i != 0 {
                write!(f, ", ")?;
            }
            if self.0[i] {
                write!(f, "1")?;
            } else {
                write!(f, "0")?;
            }
        }
        write!(f, "]")?;
        Ok(())
    }
}

impl SimpleCycle {
    pub fn new(bitset: FixedBitSet) -> SimpleCycle {
        SimpleCycle(bitset)
    }
    /// convert Cycle (a list of edges) into SimpleCycle
    pub fn from<N, E>(graph: &UnGraph<N, E>, cycle: &Cycle) -> SimpleCycle {
        let n = graph.edge_count();
        let mut bitset = FixedBitSet::with_capacity(n);

        for edge in cycle.edges().iter() {
            bitset.set(edge.index(), true);
        }

        SimpleCycle(bitset)
    }
    pub fn from_bits(n: usize, bits: &[usize]) -> SimpleCycle {
        let mut bitset = FixedBitSet::with_capacity(n);
        for &bit in bits.iter() {
            bitset.set(bit, true);
        }
        SimpleCycle(bitset)
    }
    /// access to bitset (immutable reference)
    pub fn bitset(&self) -> &FixedBitSet {
        &self.0
    }
    /// length of the bitset.
    /// Should correspond to the number of edges in the graph.
    pub fn len(&self) -> usize {
        self.0.len()
    }
    /// Get the first non-zero element in bitset.
    fn first_element(&self) -> Option<usize> {
        self.bitset().ones().next()
    }
    fn set(&mut self, index: usize, enabled: bool) {
        self.0.set(index, enabled);
    }
    ///
    /// Convert simple cycle into cycle (= list of edges)
    /// Assuming the SimpleCycle is
    /// * irreducible (the cycle is not a union of two disjoint cycle)
    /// * edge-simple (no edge is used multiple times)
    /// the traverse is easy.
    ///
    pub fn to_cycle<N, E>(&self, graph: &UnGraph<N, E>) -> Option<Cycle> {
        let mut path = Vec::new();
        let mut unvisited = self.clone();

        //            edge
        // init_node ------> node ...
        let mut edge = EdgeIndex::new(unvisited.first_element().expect("cycle is empty"));
        // update path and unvisited
        path.push(edge);
        unvisited.set(edge.index(), false);
        let (init_node, mut node) = graph.edge_endpoints(edge).unwrap();

        // pick a child of node that is in SimpleCycle
        // assert that only one child satisfies the condition.
        loop {
            let next_edges: Vec<_> = graph
                .edges(node)
                .map(|new_edge_ref| new_edge_ref.id())
                .filter(|&new_edge| new_edge != edge && unvisited.bitset()[new_edge.index()])
                .collect();

            if next_edges.len() >= 2 {
                // assert!(next_edges.len() < 2, "cycle contain a repeated node");
                return None;
            }

            // Prev iteration:
            //            edge         next_edges[0]
            // init_node ------> node --------------->
            //
            // Next iteration:
            //                             edge
            //                        ---------------> node
            if next_edges.len() == 1 {
                // len=1
                // move to next
                edge = next_edges[0];
                // update path and unvisited
                path.push(edge);
                unvisited.set(edge.index(), false);
                node = other_endpoint(graph, edge, node).expect("");
            } else {
                // len=0
                // returned to the original node
                if node != init_node {
                    panic!("could not get back to the original node");
                }
                break;
            }
        }

        Some(Cycle::new(path))
    }
    pub fn is_disjoint(&self, other: &SimpleCycle) -> bool {
        self.0.is_disjoint(&other.0)
    }
    pub fn symmetric_difference(&self, other: &SimpleCycle) -> SimpleCycle {
        let mut x = self.0.clone();
        x.symmetric_difference_with(&other.0);
        SimpleCycle::new(x)
    }
}

fn other_endpoint<N, E>(
    graph: &UnGraph<N, E>,
    edge: EdgeIndex,
    node: NodeIndex,
) -> Option<NodeIndex> {
    let (v, w) = graph.edge_endpoints(edge).unwrap();
    // if other
    if v == node && w != node {
        Some(w)
    } else if v != node && w == node {
        Some(v)
    } else if v == node && w == node {
        None
    } else {
        panic!("Neither endpoint node is not the specified node")
    }
}

impl Cycle {
    /// constructor from vec of edgeindex
    pub fn new(edges: Vec<EdgeIndex>) -> Cycle {
        Cycle(edges)
    }
    /// constructor from vec of usize slice
    pub fn from(indexes: &[usize]) -> Cycle {
        let edges = indexes.iter().map(|&i| EdgeIndex::new(i)).collect();
        Cycle::new(edges)
    }
    pub fn edges(&self) -> &[EdgeIndex] {
        &self.0
    }
    ///
    pub fn from_bitset(bitset: &FixedBitSet) -> Cycle {
        unimplemented!();
    }
    pub fn to_simple_cycle<N, E>(&self, graph: &UnGraph<N, E>) -> SimpleCycle {
        let n = graph.edge_count();
        let mut bitset = FixedBitSet::with_capacity(n);

        for edge in self.0.iter() {
            bitset.set(edge.index(), true);
        }

        SimpleCycle::new(bitset)
    }
    fn min_index(&self) -> usize {
        // find the index i (i=0,..,n-1) such that the suffix e[i:] is minimum
        let mut i0 = 0;
        let n = self.0.len();
        for i in 1..n {
            if let Ordering::Greater = cmp(&self.0, i0, i) {
                // i0 is (strictly) greater than i
                // that is the minimum candidate should be changed to i from i0
                i0 = i;
            }
        }
        i0
    }
    /// normalize the cycle
    /// so that a index vector will start in the minimum index.
    pub fn normalize(self) -> Cycle {
        let i = self.min_index();
        let mut new_cycle = self;
        new_cycle.0.rotate_left(i);
        new_cycle
    }
}

///
/// compare xs[i:] and xs[j:]
///
/// "`xs[i:]` is less/greater than `xs[j:]`?"
///
fn cmp<X: PartialOrd + Copy>(xs: &[X], i: usize, j: usize) -> Ordering {
    let n = xs.len();
    if i == j {
        return Ordering::Equal;
    }
    for k in 0..n {
        let xik = xs[(i + k) % n];
        let xjk = xs[(j + k) % n];
        if xik < xjk {
            return Ordering::Less;
        } else if xik > xjk {
            return Ordering::Greater;
        }
    }
    // all elements are the same
    Ordering::Equal
}

///
/// modify EdgeCopyNums along the cycle in +1 (if is_rev=false) or -1 (if is_rev=true)
///
/// is_rev=false the first edge will be +1.
///
pub fn apply_cycle<N, E>(
    graph: &DiGraph<N, E>,
    copy_nums: &EdgeCopyNums,
    cycle: &Cycle,
    is_rev: bool,
) -> Option<EdgeCopyNums> {
    let mut ret = copy_nums.clone();

    let mut prev_edge = None;
    let mut dir_is_rev = is_rev;

    for &edge in cycle.edges() {
        // the edge is opposite, flip the direction
        if prev_edge.is_some() && !is_same_dir(graph, prev_edge.unwrap(), edge) {
            dir_is_rev = !dir_is_rev;
        }
        if dir_is_rev {
            if ret[edge] == 0 {
                return None;
            }
            ret[edge] -= 1;
        } else {
            if ret[edge] == usize::MAX {
                return None;
            }
            ret[edge] += 1;
        }
        prev_edge = Some(edge);
    }

    Some(ret)
}

/// adjacent two edges are same direction or not?
///
/// true  if --->---> or <---<---
/// false if ---><--- or <------>
///
fn is_same_dir<N, E>(graph: &DiGraph<N, E>, e_a: EdgeIndex, e_b: EdgeIndex) -> bool {
    let (v_a, w_a) = graph.edge_endpoints(e_a).unwrap();
    let (v_b, w_b) = graph.edge_endpoints(e_b).unwrap();

    // self-loop check
    if v_a == w_a || v_b == w_b {
        panic!("direction of self-loop will be ambiguous");
    }

    if w_a == v_b || v_a == w_b {
        // --->--->
        // or
        // <---<---
        true
    } else if w_a == w_b || v_a == v_b {
        // ---><---
        // or
        // <------>
        false
    } else {
        panic!("four endpoints are not connected");
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{ei, ni};

    #[test]
    fn cycle_compare() {
        assert_eq!(cmp(&[3, 2, 4, 5, 1], 0, 1), Ordering::Greater);
        assert_eq!(cmp(&[3, 2, 4, 5, 1], 1, 0), Ordering::Less);
        assert_eq!(cmp(&[3, 2, 4, 5, 1], 0, 0), Ordering::Equal);
        assert_eq!(cmp(&[3, 2, 4, 3, 2, 4], 0, 3), Ordering::Equal);
    }

    #[test]
    fn cycle_normalize() {
        let c1 = Cycle::from(&[3, 2, 4, 5, 1]).normalize();
        println!("{:?}", c1);
        assert_eq!(c1, Cycle::from(&[1, 3, 2, 4, 5]));
        let c1 = Cycle::from(&[1, 1, 2, 1, 3]).normalize();
        println!("{:?}", c1);
        assert_eq!(c1, Cycle::from(&[1, 1, 2, 1, 3]));
        let c1 = Cycle::from(&[5, 1, 1, 2, 1, 3]).normalize();
        println!("{:?}", c1);
        assert_eq!(c1, Cycle::from(&[1, 1, 2, 1, 3, 5]));
        let c1 = Cycle::from(&[1]).normalize();
        println!("{:?}", c1);
        assert_eq!(c1, Cycle::from(&[1]));
    }

    #[test]
    fn cycle_overlap_and_unite() {
        let mut g: UnGraph<(), ()> = UnGraph::from_edges(&[
            (0, 1),
            (0, 2),
            (1, 2),
            (1, 3),
            (2, 3),
            (1, 3), // parallel edge
            (1, 1), // self loop
        ]);
        let c1 = SimpleCycle::from(&g, &Cycle::from(&[0, 1, 2]));
        let c2 = SimpleCycle::from(&g, &Cycle::from(&[2, 3, 4]));
        let c3 = SimpleCycle::from(&g, &Cycle::from(&[3, 5]));
        println!("{:?}", c1.is_disjoint(&c2));
        println!("{:?}", c1.is_disjoint(&c3));
        println!("{:?}", c2.is_disjoint(&c3));
        println!("{}", c1.symmetric_difference(&c2));

        let c1a = c1.to_cycle(&g).unwrap();
        println!("{}", c1a);
        assert_eq!(c1a, Cycle::from(&[0, 2, 1]));

        let c2a = c2.to_cycle(&g).unwrap();
        println!("{}", c2a);
        assert_eq!(c2a, Cycle::from(&[2, 4, 3]));

        let c3a = c3.to_cycle(&g).unwrap();
        println!("{}", c3a);
        assert_eq!(c3a, Cycle::from(&[3, 5]));
    }

    #[test]
    fn other_endpoint_test() {
        let g: UnGraph<(), ()> = UnGraph::from_edges(&[(0, 1), (0, 2), (1, 2), (3, 3)]);
        assert_eq!(other_endpoint(&g, ei(0), ni(0)), Some(ni(1)));
        assert_eq!(other_endpoint(&g, ei(0), ni(1)), Some(ni(0)));
        assert_eq!(other_endpoint(&g, ei(1), ni(0)), Some(ni(2)));
        assert_eq!(other_endpoint(&g, ei(1), ni(2)), Some(ni(0)));
        assert_eq!(other_endpoint(&g, ei(1), ni(2)), Some(ni(0)));
        assert_eq!(other_endpoint(&g, ei(3), ni(3)), None);
    }
    #[test]
    #[should_panic]
    fn other_endpoint_test_panic() {
        let g: UnGraph<(), ()> = UnGraph::from_edges(&[(0, 1), (0, 2), (1, 2), (3, 3)]);
        other_endpoint(&g, ei(0), ni(3));
    }

    #[test]
    fn apply_cycle_test() {
        //
        // 0 -> 1 -> 2 -> 3
        //  \           /
        //   -> 4 -> 5 >
        //
        let g: DiGraph<(), ()> =
            DiGraph::from_edges(&[(0, 1), (1, 2), (2, 3), (0, 4), (4, 5), (5, 3)]);
        let n = EdgeCopyNums::from_slice(&[1, 1, 1, 1, 1, 1], 0);
        let c = Cycle::from(&[0, 1, 2, 5, 4, 3]);
        println!("n={}", n);
        let n1 = apply_cycle(&g, &n, &c, false).unwrap();
        println!("n1={}", n1);
        assert_eq!(n1.to_vec(), vec![2, 2, 2, 0, 0, 0]);
        let n2 = apply_cycle(&g, &n, &c, true).unwrap();
        println!("n2={}", n2);
        assert_eq!(n2.to_vec(), vec![0, 0, 0, 2, 2, 2]);
        let n3 = apply_cycle(&g, &n2, &c, false).unwrap();
        assert_eq!(n3.to_vec(), vec![1, 1, 1, 1, 1, 1]);
        let n4 = apply_cycle(&g, &n2, &c, true);
        assert!(n4.is_none());
    }
}
