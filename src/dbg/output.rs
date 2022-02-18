//!
//! Output related functions of Dbg
//!
use super::dbg::{Dbg, DbgEdge, DbgNode};
use crate::io::cytoscape::{EdgeAttr, ElementV2, NodeAttr};
use petgraph::dot::Dot;

impl<N, E> std::fmt::Display for Dbg<N, E>
where
    N: DbgNode + std::fmt::Display,
    E: DbgEdge + std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", Dot::with_config(&self.graph, &[]))
    }
}

//
// cytoscape output
//
impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    /// Convert de bruijn graph into
    /// collection of elements for cytoscape parsing
    fn to_cytoscape_elements(&self) -> Vec<ElementV2<N::Kmer>> {
        let mut elements = Vec::new();

        // add nodes
        for (node, weight) in self.nodes() {
            elements.push(ElementV2::Node {
                id: node,
                label: Some(weight.kmer().clone()),
                attrs: vec![NodeAttr::CopyNum(weight.copy_num())],
            });
        }

        // add edges
        let n = self.n_nodes();
        for (edge, s, t, weight) in self.edges() {
            elements.push(ElementV2::Edge {
                // element id should be unique, for all nodes and edges
                id: edge,
                source: s,
                target: t,
                label: None,
                attrs: vec![],
            });
        }
        elements
    }
    /// Convert de bruijn graph into
    /// cytoscape-parsable JSON string
    pub fn to_cytoscape(&self) -> String {
        let elements = self.to_cytoscape_elements();
        serde_json::to_string(&elements).unwrap()
    }
    // /// output as cytoscape json format
    // /// by converting into edge centric dbg
    // pub fn to_edge_cytoscape(&self) -> String {
    //     self.to_edbg().to_cytoscape()
    // }
}
