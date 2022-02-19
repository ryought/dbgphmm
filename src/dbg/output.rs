//!
//! Output related functions of Dbg
//!
use super::dbg::{Dbg, DbgEdge, DbgNode};
use crate::common::{ei, ni};
use crate::io::cytoscape::{EdgeAttr, EdgeAttrVec, ElementV2, NodeAttr, NodeAttrVec};
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
    fn to_cytoscape_elements_with_attrs(
        &self,
        node_attrs: &[NodeAttrVec],
        edge_attrs: &[EdgeAttrVec],
    ) -> Vec<ElementV2<N::Kmer>> {
        let mut elements = Vec::new();

        // add nodes
        for (node, weight) in self.nodes() {
            // attrs for the node
            let mut attrs = vec![NodeAttr::CopyNum(weight.copy_num())];
            for node_attr in node_attrs {
                match node_attr {
                    NodeAttrVec::CopyNum(v) => attrs.push(NodeAttr::CopyNum(v[node])),
                    NodeAttrVec::Freq(v) => attrs.push(NodeAttr::Freq(v[node])),
                }
            }
            elements.push(ElementV2::Node {
                id: node,
                label: Some(weight.kmer().clone()),
                attrs,
            });
        }

        // add edges
        let n = self.n_nodes();
        for (edge, s, t, weight) in self.edges() {
            // attrs for the edge
            let mut attrs = vec![];
            for edge_attr in edge_attrs {
                match edge_attr {
                    EdgeAttrVec::TrueCopyNum(v) => attrs.push(EdgeAttr::TrueCopyNum(v[edge])),
                    EdgeAttrVec::Freq(v) => attrs.push(EdgeAttr::Freq(v[edge])),
                }
            }
            elements.push(ElementV2::Edge {
                // element id should be unique, for all nodes and edges
                id: edge,
                source: s,
                target: t,
                label: None,
                attrs,
            });
        }
        elements
    }
    /// Convert de bruijn graph into
    /// collection of elements for cytoscape parsing
    /// without any node/edge attributes
    fn to_cytoscape_elements(&self) -> Vec<ElementV2<N::Kmer>> {
        self.to_cytoscape_elements_with_attrs(&[], &[])
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dbg::mocks::mock_simple;
    use crate::vector::{DenseStorage, EdgeVec, NodeVec};

    #[test]
    fn dbg_mock_simple_cytoscape() {
        let dbg = mock_simple();
        // (1) with no attrs
        let json = dbg.to_cytoscape();
        println!("{}", json);
        let json_true = r#"[{"group":"nodes","data":{"id":"n0","label":"TTNN","attrs":[{"type":"copy_num","value":1}]}},{"group":"nodes","data":{"id":"n1","label":"TNNN","attrs":[{"type":"copy_num","value":1}]}},{"group":"nodes","data":{"id":"n2","label":"AGCT","attrs":[{"type":"copy_num","value":1}]}},{"group":"nodes","data":{"id":"n3","label":"GCTT","attrs":[{"type":"copy_num","value":1}]}},{"group":"nodes","data":{"id":"n4","label":"AAAG","attrs":[{"type":"copy_num","value":1}]}},{"group":"nodes","data":{"id":"n5","label":"CTTG","attrs":[{"type":"copy_num","value":1}]}},{"group":"nodes","data":{"id":"n6","label":"ATTN","attrs":[{"type":"copy_num","value":1}]}},{"group":"nodes","data":{"id":"n7","label":"TTGA","attrs":[{"type":"copy_num","value":1}]}},{"group":"nodes","data":{"id":"n8","label":"NAAA","attrs":[{"type":"copy_num","value":1}]}},{"group":"nodes","data":{"id":"n9","label":"NNAA","attrs":[{"type":"copy_num","value":1}]}},{"group":"nodes","data":{"id":"n10","label":"NNNA","attrs":[{"type":"copy_num","value":1}]}},{"group":"nodes","data":{"id":"n11","label":"AAGC","attrs":[{"type":"copy_num","value":1}]}},{"group":"nodes","data":{"id":"n12","label":"TGAT","attrs":[{"type":"copy_num","value":1}]}},{"group":"nodes","data":{"id":"n13","label":"GATT","attrs":[{"type":"copy_num","value":1}]}},{"group":"edges","data":{"id":"e0","source":"n0","target":"n1","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e1","source":"n1","target":"n10","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e2","source":"n2","target":"n3","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e3","source":"n3","target":"n5","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e4","source":"n4","target":"n11","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e5","source":"n5","target":"n7","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e6","source":"n6","target":"n0","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e7","source":"n7","target":"n12","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e8","source":"n8","target":"n4","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e9","source":"n9","target":"n8","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e10","source":"n10","target":"n9","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e11","source":"n11","target":"n2","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e12","source":"n12","target":"n13","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e13","source":"n13","target":"n6","label":null,"attrs":[]}}]"#;
        assert_eq!(json, json_true);

        // (2) with attrs
        let mut v1 = NodeVec::new(dbg.n_nodes(), 1);
        v1[ni(0)] = 100;
        let n1 = NodeAttrVec::CopyNum(v1);
        let mut v2 = NodeVec::new(dbg.n_nodes(), 0.0);
        v2[ni(1)] = 10.1;
        let n2 = NodeAttrVec::Freq(v2);
        let mut w1 = EdgeVec::new(dbg.n_nodes(), 1);
        w1[ei(0)] = 100;
        let e1 = EdgeAttrVec::TrueCopyNum(w1);
        let elements = dbg.to_cytoscape_elements_with_attrs(&[n1, n2], &[e1]);
        // let json = serde_json::to_string_pretty(&elements).unwrap();
        let json = serde_json::to_string(&elements).unwrap();
        println!("{}", json);
        let json_true = r#"[{"group":"nodes","data":{"id":"n0","label":"TTNN","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":100},{"type":"freq","value":0.0}]}},{"group":"nodes","data":{"id":"n1","label":"TNNN","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":10.1}]}},{"group":"nodes","data":{"id":"n2","label":"AGCT","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}]}},{"group":"nodes","data":{"id":"n3","label":"GCTT","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}]}},{"group":"nodes","data":{"id":"n4","label":"AAAG","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}]}},{"group":"nodes","data":{"id":"n5","label":"CTTG","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}]}},{"group":"nodes","data":{"id":"n6","label":"ATTN","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}]}},{"group":"nodes","data":{"id":"n7","label":"TTGA","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}]}},{"group":"nodes","data":{"id":"n8","label":"NAAA","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}]}},{"group":"nodes","data":{"id":"n9","label":"NNAA","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}]}},{"group":"nodes","data":{"id":"n10","label":"NNNA","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}]}},{"group":"nodes","data":{"id":"n11","label":"AAGC","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}]}},{"group":"nodes","data":{"id":"n12","label":"TGAT","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}]}},{"group":"nodes","data":{"id":"n13","label":"GATT","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}]}},{"group":"edges","data":{"id":"e0","source":"n0","target":"n1","label":null,"attrs":[{"type":"true_copy_num","value":100}]}},{"group":"edges","data":{"id":"e1","source":"n1","target":"n10","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e2","source":"n2","target":"n3","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e3","source":"n3","target":"n5","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e4","source":"n4","target":"n11","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e5","source":"n5","target":"n7","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e6","source":"n6","target":"n0","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e7","source":"n7","target":"n12","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e8","source":"n8","target":"n4","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e9","source":"n9","target":"n8","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e10","source":"n10","target":"n9","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e11","source":"n11","target":"n2","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e12","source":"n12","target":"n13","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e13","source":"n13","target":"n6","label":null,"attrs":[{"type":"true_copy_num","value":1}]}}]"#;
        assert_eq!(json, json_true);
    }
}
