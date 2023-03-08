//!
//! Cycle in graph
//!
use crate::dbg::dbg::{EdgeCopyNums, NodeCopyNums};
pub use crate::graph_public::cycle::Cycle;
use crate::utils::breakpoints;
use fixedbitset::FixedBitSet;
use itertools::Itertools;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex, UnGraph};
use petgraph::visit::EdgeRef;

///
/// cycle with direction info (as a list of (edge, direction))
///
#[derive(Debug, Clone, PartialEq)]
pub struct CycleWithDir(Vec<(EdgeIndex, bool)>);
impl CycleWithDir {
    pub fn empty() -> Self {
        Self(Vec::new())
    }
    pub fn new(edges: Vec<(EdgeIndex, bool)>) -> Self {
        Self(edges)
    }
    pub fn edges(&self) -> &[(EdgeIndex, bool)] {
        &self.0
    }
    /// reverse all directions in the cycle
    pub fn reverse_dir(self) -> Self {
        Self(self.0.into_iter().map(|(e, is_rev)| (e, !is_rev)).collect())
    }
    ///
    /// List up breakpoints of the cycle.
    ///
    /// breakpoint is index `i` such that i-th edge and i+1-th edge have different direction.
    ///
    fn breakpoints(&self) -> Vec<usize> {
        breakpoints(&self.0)
    }
    ///
    /// collapse the adj edges with same direction
    /// `[e1+,e2+,e3+,e4-,e5-,e6-] -> [([e1,e2,e3],+),([e4,e5,e6],-)]`
    ///
    pub fn collapse_dir(&self) -> Vec<(Vec<EdgeIndex>, bool)> {
        if self.edges().len() == 0 {
            // empty cycle
            Vec::new()
        } else {
            let breakpoints = self.breakpoints();
            let n = breakpoints.len();
            let mut ret = Vec::new();
            // is the same dir only
            // look for the breaking point
            if n == 0 {
                ret.push(to_collapsed(&self.0));
            } else {
                for i in 0..n {
                    let breakpoint_a = breakpoints[i];
                    let breakpoint_b = breakpoints[(i + 1) % n];
                    let segment = get_round_slice(&self.0, breakpoint_a, breakpoint_b);
                    ret.push(to_collapsed(&segment));
                }
            }
            ret
        }
    }
}
///
/// get `xs[i..j]` with rounding index `i,j (0 <= i,j < len(xs))`
///
/// * i <= j then usual xs[i..j]
/// * i > j  then xs[i..] xs[..j] is returned
///
fn get_round_slice<X: Clone>(xs: &[X], i: usize, j: usize) -> Vec<X> {
    let n = xs.len();
    assert!(i < n);
    assert!(j < n);
    if i <= j {
        xs[i..j].to_vec()
    } else {
        [&xs[i..], &xs[..j]].concat()
    }
}
fn get_all_same_value<T: PartialEq + Clone>(slice: &[T]) -> Option<T> {
    match slice.get(0) {
        Some(first) => {
            if slice.iter().all(|x| x == first) {
                Some(first.clone())
            } else {
                None
            }
        }
        None => None,
    }
}
///
/// Vec<(T, bool)> into (Vec<T>, bool) assume bool has the same entry.
///
fn to_collapsed<T: Clone>(slice: &[(T, bool)]) -> (Vec<T>, bool) {
    let xs: Vec<T> = slice.iter().map(|(x, _)| x.clone()).collect();
    let bs: Vec<bool> = slice.iter().map(|(_, b)| *b).collect();
    let b = get_all_same_value(&bs).unwrap();
    (xs, b)
}
impl std::fmt::Display for CycleWithDir {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.0
                .iter()
                .map(|(e, is_rev)| format!("{}{}", e.index(), if *is_rev { "-" } else { "+" }))
                .join(",")
        )
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

pub fn apply_cycle_with_dir(
    copy_nums: &EdgeCopyNums,
    cycle: &CycleWithDir,
) -> Option<EdgeCopyNums> {
    let mut ret = copy_nums.clone();
    for &(edge, dir_is_rev) in cycle.edges() {
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
    }
    Some(ret)
}

pub fn to_cycle_with_dir<N, E>(graph: &DiGraph<N, E>, cycle: &Cycle) -> CycleWithDir {
    let mut prev_edge = None;
    let mut dir_is_rev = false; // direction of the first edge is always forward.
    let mut ret = Vec::new();

    for &edge in cycle.edges() {
        // the edge is opposite, flip the direction
        if prev_edge.is_some() && !is_same_dir(graph, prev_edge.unwrap(), edge) {
            dir_is_rev = !dir_is_rev;
        }
        ret.push((edge, dir_is_rev));
        prev_edge = Some(edge);
    }

    CycleWithDir::new(ret)
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
    let cycle_with_dir = to_cycle_with_dir(graph, cycle);
    if is_rev {
        apply_cycle_with_dir(copy_nums, &cycle_with_dir.reverse_dir())
    } else {
        apply_cycle_with_dir(copy_nums, &cycle_with_dir)
    }
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

    #[test]
    fn get_round_slice_test() {
        let xs = vec![0, 1, 2, 3, 4, 5];
        println!("{:?}", get_round_slice(&xs, 4, 1));
        assert_eq!(get_round_slice(&xs, 0, 0).len(), 0);
        assert_eq!(get_round_slice(&xs, 2, 2).len(), 0);
        assert_eq!(get_round_slice(&xs, 1, 4), vec![1, 2, 3]);
        assert_eq!(get_round_slice(&xs, 1, 5), vec![1, 2, 3, 4]);
        assert_eq!(get_round_slice(&xs, 4, 1), vec![4, 5, 0]);
    }
    #[test]
    fn get_all_same_value_test() {
        assert_eq!(get_all_same_value(&[0, 0, 0, 0]), Some(0));
        assert_eq!(get_all_same_value(&[0, 1, 0, 0]), None);
        assert_eq!(get_all_same_value::<usize>(&[]), None);
    }

    #[test]
    fn cycle_with_dir_test() {
        {
            let c0 = CycleWithDir::new(vec![
                (ei(0), false),
                (ei(1), true),
                (ei(2), true),
                (ei(3), true),
                (ei(4), false),
                (ei(5), false),
            ]);
            // c0 has two breakpoints
            // - between c0[0] and c0[1]
            // - between c0[3] and c0[4]
            println!("{:?}", c0.breakpoints());
            println!("{:?}", c0.collapse_dir());
            assert_eq!(c0.breakpoints(), vec![1, 4]);
            assert_eq!(
                c0.collapse_dir(),
                vec![
                    (vec![ei(1), ei(2), ei(3)], true),
                    (vec![ei(4), ei(5), ei(0)], false),
                ]
            );
        }

        {
            let c0 = CycleWithDir::new(vec![
                (ei(0), true),
                (ei(1), true),
                (ei(2), true),
                (ei(3), false),
                (ei(4), false),
                (ei(5), false),
            ]);
            assert_eq!(c0.breakpoints(), vec![0, 3]);
            assert_eq!(
                c0.collapse_dir(),
                vec![
                    (vec![ei(0), ei(1), ei(2)], true),
                    (vec![ei(3), ei(4), ei(5)], false),
                ]
            );
        }

        {
            let c0 = CycleWithDir::new(vec![
                (ei(0), true),
                (ei(1), true),
                (ei(2), true),
                (ei(3), true),
            ]);
            assert_eq!(c0.breakpoints().len(), 0);
            assert_eq!(
                c0.collapse_dir(),
                vec![(vec![ei(0), ei(1), ei(2), ei(3)], true)]
            );
        }

        {
            let c0 = CycleWithDir::new(vec![
                (ei(0), true),
                (ei(1), true),
                (ei(2), false),
                (ei(3), false),
                (ei(4), true),
                (ei(5), true),
                (ei(6), false),
                (ei(7), false),
                (ei(8), true),
            ]);
            assert_eq!(c0.breakpoints(), vec![2, 4, 6, 8]);
            assert_eq!(
                c0.collapse_dir(),
                vec![
                    (vec![ei(2), ei(3)], false),
                    (vec![ei(4), ei(5)], true),
                    (vec![ei(6), ei(7)], false),
                    (vec![ei(8), ei(0), ei(1)], true),
                ]
            );
        }
    }
}
