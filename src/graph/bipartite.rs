//!
//! Constructor
//!
//! * use function
//! * iterators on in_nodes out_nodes edges
//! * access by index
//!

///
/// Struct representing a bipartite graph
///
/// ```text
/// in_nodes:  AXXX, CXXX, GXXX, TXXX
/// out_nodes: XXXA, XXXC, XXXG, XXXT
/// ```
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
        edges: F,
    ) -> Self {
        unimplemented!();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn starting() {
        assert!();
        assert_eq!();
        assert_ne!();
    }
}
