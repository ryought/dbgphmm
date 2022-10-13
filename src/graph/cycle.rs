//!
//! Cycle in graph
//!
use fixedbitset::FixedBitSet;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex, UnGraph};
use std::cmp::Ordering;

///
/// cycle (as a list of edges)
///
#[derive(Debug, Clone, PartialEq)]
pub struct Cycle(Vec<EdgeIndex>);

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
    pub fn to_cycle<N, E>(graph: &UnGraph<N, E>) -> Cycle {
        unimplemented!();
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

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

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
    }
}
