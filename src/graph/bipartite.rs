//!
//! `Bipartite`
//!
//! # method requirements
//!
//! * use function
//! * iterators on in_nodes out_nodes edges
//! * access by index
//!
// use arrayvec::ArrayVec;
use itertools::iproduct;

///
/// Struct representing a bipartite graph
///
/// ```text
/// in_nodes:  AXXX, CXXX, GXXX, TXXX
/// out_nodes: XXXA, XXXC, XXXG, XXXT
/// ```
///
/// ## TODO
///
/// * change Vec into ArrayVec to avoid heap allocation.
///
#[derive(Debug, Clone)]
pub struct Bipartite<I, N, E> {
    pub id: I,
    pub in_nodes: Vec<N>,
    pub out_nodes: Vec<N>,
    pub edges: Vec<E>,
}

impl<I, N, E> Bipartite<I, N, E> {
    pub fn new(id: I, in_nodes: Vec<N>, out_nodes: Vec<N>, edges: Vec<E>) -> Self {
        Bipartite {
            id,
            in_nodes,
            out_nodes,
            edges,
        }
    }
    pub fn from<F: Fn(usize, usize) -> E>(
        id: I,
        in_nodes: Vec<N>,
        out_nodes: Vec<N>,
        edge_fn: F,
    ) -> Self {
        let edges = iproduct!(0..in_nodes.len(), 0..out_nodes.len())
            .map(|(i, j)| edge_fn(i, j))
            .collect();
        Bipartite::new(id, in_nodes, out_nodes, edges)
    }
    pub fn n_in(&self) -> usize {
        self.in_nodes.len()
    }
    pub fn n_out(&self) -> usize {
        self.out_nodes.len()
    }
    ///
    /// Get an in-node by specifying index.
    ///
    pub fn in_node(&self, index: usize) -> &N {
        assert!(index < self.n_in());
        &self.in_nodes[index]
    }
    ///
    /// Get an out-node by specifying index.
    ///
    pub fn out_node(&self, index: usize) -> &N {
        assert!(index < self.n_out());
        &self.out_nodes[index]
    }
    ///
    /// Get an edge between in_node and out_node, specified by two indexs.
    ///
    /// ## Implementation details
    ///
    /// For example, in the case of `n_in=2, n_out=3`,
    /// the edges are stored as `[a00, a01, a02, a10, a11, a12]`
    /// (first index is `in` and the last index is `out`).
    ///
    pub fn edge(&self, index_in: usize, index_out: usize) -> &E {
        assert!(index_in < self.n_in());
        assert!(index_out < self.n_out());
        &self.edges[index_in * self.n_out() + index_out]
    }
    pub fn iter_in_nodes(&self) -> impl Iterator<Item = &N> + '_ {
        self.in_nodes.iter()
    }
    pub fn iter_out_nodes(&self) -> impl Iterator<Item = &N> + '_ {
        self.out_nodes.iter()
    }
}

impl<I, N, E> std::fmt::Display for Bipartite<I, N, E>
where
    I: std::fmt::Display,
    N: std::fmt::Display,
    E: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "id\t{}\n", self.id)?;
        for i in 0..self.n_in() {
            write!(f, "in\t{}\t{}\n", i, self.in_node(i))?;
        }
        for j in 0..self.n_out() {
            write!(f, "out\t{}\t{}\n", j, self.out_node(j))?;
        }
        for (i, j) in iproduct!(0..self.n_in(), 0..self.n_out()) {
            write!(f, "in{}\tout{}\t{}\n", i, j, self.edge(i, j))?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bipartite_construct() {
        let in_nodes = vec![1, 2, 3, 4];
        let out_nodes = vec![11, 12, 13];
        let bi: Bipartite<usize, usize, usize> =
            Bipartite::from(111, in_nodes.clone(), out_nodes.clone(), |i, j| {
                in_nodes[i] + out_nodes[j]
            });
        assert_eq!(bi.id, 111);
        assert_eq!(bi.in_nodes, vec![1, 2, 3, 4]);
        assert_eq!(bi.out_nodes, vec![11, 12, 13]);
        assert_eq!(
            bi.edges,
            vec![12, 13, 14, 13, 14, 15, 14, 15, 16, 15, 16, 17]
        );

        assert_eq!(*bi.in_node(0), 1);
        assert_eq!(*bi.in_node(3), 4);
        assert_eq!(*bi.out_node(0), 11);
        assert_eq!(*bi.out_node(2), 13);
        assert_eq!(*bi.edge(0, 1), 13);
        assert_eq!(*bi.edge(3, 2), 17);

        println!("{}", bi);
    }
    #[test]
    fn iproduct_test() {
        let prod: Vec<(usize, usize)> = iproduct!(0..2, 0..3).collect();
        println!("{:?}", prod);
        assert_eq!(prod, vec![(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]);
    }
}
