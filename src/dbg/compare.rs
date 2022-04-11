//!
//! Comparison methods of two de Bruijn graphs.
//!
use super::dbg::{Dbg, DbgEdge, DbgNode};

impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {}

//
// tests
//
#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::sequence_to_string;

    #[test]
    fn dbg_compare_simple() {
        /*
        let dbg = mock_simple();
        let seqs = dbg.to_seqs();
        println!("{}", dbg);
        for seq in seqs.iter() {
            println!("{}", sequence_to_string(seq));
        }
        assert_eq!(seqs, vec![b"NNAAAGCTTGATTN"]);
        */
    }
}
