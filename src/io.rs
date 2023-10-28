pub mod json;

use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use petgraph::graph::{DefaultIx, DiGraph, EdgeIndex, IndexType, NodeIndex};
use std::fs::File;
use std::io::prelude::*;

///
/// write string into a file
///
/// ```text
/// write_string("hoge.txt", "hogehogehoge\nhogehoge");
/// ```
///
pub fn write_string(filename: &str, string: &str) -> std::io::Result<()> {
    let mut file = File::create(&filename)?;
    file.write_all(string.as_bytes())?;
    Ok(())
}

#[derive(Clone, PartialEq, PartialOrd, Eq, Copy, Default, Hash, Ord, Debug)]
pub struct CompactIndex(usize);

unsafe impl IndexType for CompactIndex {
    #[inline(always)]
    fn new(x: usize) -> Self {
        CompactIndex(x)
    }
    #[inline(always)]
    fn index(&self) -> usize {
        self.0
    }
    #[inline(always)]
    fn max() -> Self {
        CompactIndex(std::usize::MAX)
    }
}

pub fn f() -> DiGraph<usize, usize, CompactIndex> {
    let mut g: DiGraph<usize, usize, CompactIndex> = DiGraph::default();
    let v = g.add_node(10);
    let w = g.add_node(2);
    g.add_edge(v, w, 100);

    for e in g.edge_indices() {
        println!("e={:?}", e);
    }
    g
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gz_test() {
        f();
        /*
        let vec = Vec::new();
        let mut e = GzEncoder::new(vec, Compression::default());
        e.write_all(b"hoge");
        e.finish().unwrap();

        println!("vec={:?}", vec);

        let mut d = GzDecoder::new(vec);
        let mut s = String::new();
        d.read_to_string(&mut s).unwrap();
        println!("{}", s);
        */
    }
}
