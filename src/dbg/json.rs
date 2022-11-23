//!
//! Dbg serialization
//!
use super::dbg::{Dbg, DbgEdge, DbgNode};
use crate::common::{CopyNum, Seq, Sequence};
use crate::kmer::{KmerLike, NullableKmer};
use petgraph::graph::{DiGraph, NodeIndex};
use serde::{Deserialize, Serialize};
use serde_with::serde_as;
use std::path::Path;
use std::{fmt, fs, io};

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DbgAsJson {
    k: usize,
    genome_size: usize,
    n_nodes: usize,
    n_edges: usize,
    nodes: Vec<(usize, String, CopyNum)>,
    edges: Vec<(usize, usize, usize)>,
}

impl<N, E> Dbg<N, E>
where
    N: DbgNode,
    E: DbgEdge,
{
    ///
    ///
    ///
    pub fn to_json_struct(&self) -> DbgAsJson {
        let nodes = self
            .nodes()
            .map(|(node, node_weight)| {
                (
                    node.index(),
                    node_weight.kmer().to_bases().to_str().to_owned(),
                    node_weight.copy_num(),
                )
            })
            .collect();
        let edges = self
            .edges()
            .map(|(edge, source, target, _)| (edge.index(), source.index(), target.index()))
            .collect();

        DbgAsJson {
            k: self.k(),
            genome_size: self.genome_size(),
            n_nodes: self.n_nodes(),
            n_edges: self.n_edges(),
            nodes,
            edges,
        }
    }
    ///
    ///
    ///
    pub fn from_json_struct(j: DbgAsJson) -> Self {
        let mut graph = DiGraph::new();
        for (id, kmer_string, copy_num) in j.nodes {
            let kmer = N::Kmer::from_str(&kmer_string);
            let node = graph.add_node(N::new(kmer, copy_num));
            assert_eq!(node.index(), id);
        }
        for (id, source, target) in j.edges {
            let edge = graph.add_edge(NodeIndex::new(source), NodeIndex::new(target), E::new(None));
            assert_eq!(edge.index(), id);
        }

        Dbg::from_digraph(j.k, graph)
    }
    pub fn to_json(&self) -> String {
        serde_json::to_string(&self.to_json_struct()).unwrap()
    }
    // file related
    pub fn to_json_file<P: AsRef<Path>>(&self, path: P) {
        let mut file = fs::File::create(path).unwrap();
        serde_json::to_writer(&mut file, &self.to_json_struct()).unwrap()
    }
    pub fn from_json_file<P: AsRef<Path>>(path: P) -> Self {
        let mut file = fs::File::open(path).unwrap();
        let json_struct: DbgAsJson = serde_json::from_reader(&mut file).unwrap();
        Self::from_json_struct(json_struct)
    }
}

impl<N, E> Dbg<N, E>
where
    N: DbgNode,
    E: DbgEdge,
{
    pub fn to_gfa<W: fmt::Write>(&self, f: &mut W) -> fmt::Result {
        writeln!(f, "K\t{}", self.k())?;
        writeln!(f, "G\t{}", self.genome_size())?;
        writeln!(f, "N\t{}", self.n_nodes())?;
        writeln!(f, "M\t{}", self.n_edges())?;
        for (node, node_weight) in self.nodes() {
            writeln!(
                f,
                "N\t{}\t{}\t{}",
                node.index(),
                node_weight.kmer(),
                node_weight.copy_num()
            )?;
        }
        for (edge, source, target, weight) in self.edges() {
            writeln!(
                f,
                "E\t{}\t{}\t{}\t{:?}",
                edge.index(),
                source.index(),
                target.index(),
                weight.copy_num(),
            )?;
        }
        Ok(())
    }
    pub fn to_gfa_string(&self) -> String {
        let mut output = String::new();
        self.to_gfa(&mut output).unwrap();
        output
    }
    pub fn to_gfa_file(&self, filename: &str) -> io::Result<()> {
        let s = self.to_gfa_string();
        fs::write(filename, &s)
    }
    pub fn from_gfa_str(s: &str) -> Self {
        let mut k = 0;
        for line in s.lines() {
            if line.is_empty() {
                continue;
            }
            match line.chars().nth(0).unwrap() {
                'K' => {}
                _ => panic!(),
            }
        }
        unimplemented!();
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dbg::impls::{SimpleDbg, SimpleDbgEdge, SimpleDbgNode};
    use crate::dbg::mocks;
    use crate::kmer::veckmer::{kmer, VecKmer};
    #[test]
    fn dbg_gfa_intersection_small() {
        let dbg = mocks::mock_intersection_small();
        let gfa = dbg.to_gfa_string();
        println!("{}", gfa);
        // let dbg2: SimpleDbg<VecKmer> = SimpleDbg::from_gfa_str(&gfa);
    }
    #[test]
    fn serde_json_flatten_test() {}
    #[test]
    fn dbg_json_intersection_small() {
        let dbg = mocks::mock_intersection_small();
        let nodes: Vec<_> = dbg.nodes().collect();
        let edges: Vec<_> = dbg.edges().collect();
        let dbg_as_json = dbg.to_json_struct();
        let json = serde_json::to_string(&dbg_as_json).unwrap();
        println!("{}", json);

        let d2: DbgAsJson = serde_json::from_str(&json).unwrap();
        let dbg2: SimpleDbg<VecKmer> = SimpleDbg::from_json_struct(d2);
        let nodes2: Vec<_> = dbg2.nodes().collect();
        let edges2: Vec<_> = dbg2.edges().collect();
        assert_eq!(nodes, nodes2);
        assert_eq!(edges, edges2);
    }
}
