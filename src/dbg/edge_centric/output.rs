use super::{EDbg, EDbgEdge, EDbgNode};
use crate::io::cytoscape::Element;
use crate::kmer::{KmerLike, VecKmer};
use petgraph::dot::Dot;

///
/// Normal DOT output
///
impl<N, E> std::fmt::Display for EDbg<N, E>
where
    N: EDbgNode + std::fmt::Display,
    E: EDbgEdge + std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", Dot::with_config(&self.graph, &[]))
    }
}

//
// cytoscape output
//
impl<N: EDbgNode, E: EDbgEdge> EDbg<N, E> {
    /// Convert edge-centric de bruijn graph into
    /// collection of elements for cytoscape parsing
    fn to_cytoscape_elements(&self) -> Vec<Element> {
        let mut elements = Vec::new();

        // add nodes
        for (node, weight) in self.nodes() {
            elements.push(Element::Node {
                id: node.index(),
                label: VecKmer::from_bases(b""),
            });
        }

        // add edges
        let n = self.n_nodes();
        for (edge, s, t, weight) in self.edges() {
            elements.push(Element::Edge {
                // element id should be unique, for all nodes and edges
                id: n + edge.index(),
                source: s.index(),
                target: t.index(),
                label: VecKmer::from_bases(&weight.kmer().to_bases()),
                widths: vec![],
                true_width: Some(weight.copy_num() as u32),
            });
        }
        elements
    }
    /// Convert edge-centric de bruijn graph into
    /// cytoscape-parsable JSON string
    pub fn to_cytoscape(&self) -> String {
        let elements = self.to_cytoscape_elements();
        serde_json::to_string(&elements).unwrap()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dbg::mocks::mock_simple;

    #[test]
    fn edbg_cytoscape() {
        let dbg = mock_simple();
        let json = dbg.to_edbg().to_cytoscape();
        println!("{}", json);
    }
}
