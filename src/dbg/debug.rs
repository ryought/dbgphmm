//!
//! Dbg debug visualizers
//!
use super::dbg::{Dbg, DbgEdge, DbgNode};
use crate::vector::{EdgeVec, NodeVec, Storage};
use itertools::Itertools;

impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    ///
    /// visualize node vec
    ///
    pub fn draw_with_vecs<S>(&self, nvs: &[&NodeVec<S>], evs: &[&EdgeVec<S>])
    where
        S: Storage,
        S::Item: std::fmt::Display,
    {
        for (path_id, path) in self.traverse_all().enumerate() {
            println!("#### path{} ####", path_id);
            let n = path.len();
            for i in 0..n {
                let node = path[i];
                let next_node = path[(i + 1) % n];
                let edge = self.find_edge(node, next_node).unwrap();

                let node_weight = self.node(node);
                let edge_weight = self.edge(edge);

                // show node infos
                let node_infos = nvs.iter().map(|nv| nv[node].to_string()).join("\t");
                println!(
                    "v{}\t{}\tx{}\t>{}\t{}",
                    node.index(),
                    node_weight.kmer(),
                    node_weight.copy_num(),
                    self.childs(node).count(),
                    node_infos
                );

                // show edge infos
                let edge_infos = evs.iter().map(|ev| ev[edge].to_string()).join("\t");
                println!(
                    "e{}\t|\tx{:?}\t\t{}",
                    edge.index(),
                    edge_weight.copy_num(),
                    edge_infos
                );
            }
        }
    }
}
