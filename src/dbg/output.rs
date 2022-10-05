//!
//! Output related functions of Dbg
//!
use super::dbg::{Dbg, DbgEdge, DbgNode};
use super::hashdbg_v2::HashDbg;
use crate::common::{ei, ni, StyledSequence, StyledSequenceParseError};
use crate::io::cytoscape::{EdgeAttr, EdgeAttrVec, ElementV2, NodeAttr, NodeAttrVec};
use crate::vector::{DenseStorage, EdgeVec, NodeVec};
use itertools::Itertools;
use petgraph::dot::Dot;
use std::str::FromStr;

//
// FromStr/Display implementations
//

impl<N, E> FromStr for Dbg<N, E>
where
    N: DbgNode,
    E: DbgEdge,
{
    type Err = StyledSequenceParseError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts: Vec<&str> = s.split(',').collect();

        // parse k-mer size
        let k = parts[0]
            .parse::<usize>()
            .expect("serialized Dbg is invalid in k-mer size");

        // create hashdbg and insert all reads
        let mut hd = HashDbg::new(k);
        for &part in parts[1..].iter() {
            let seq = StyledSequence::from_str(part)
                .expect("serialized Dbg is invalid in styled sequence format");
            hd.add_styled_sequence(&seq);
        }

        Ok(Self::from_hashdbg(&hd))
    }
}

///
/// Serialize Dbg as `k,(eulerian traversed sequences with style, separated with ',')`
///
/// Example:
/// ```text
/// 8,L:ATCGATCG,L:ATTTAC
/// ```
///
impl<N, E> std::fmt::Display for Dbg<N, E>
where
    N: DbgNode + std::fmt::Display,
    E: DbgEdge + std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let seqs = self
            .to_styled_seqs()
            .iter()
            .map(|seq| seq.to_string())
            .join(",");
        write!(f, "{},{}", self.k(), seqs)
    }
}

// debug output
impl<N, E> Dbg<N, E>
where
    N: DbgNode + std::fmt::Display,
    E: DbgEdge + std::fmt::Display,
{
    pub fn to_dot(&self) -> String {
        format!("{}", Dot::with_config(&self.graph, &[]))
    }
}

//
// cytoscape output
//
impl<N: DbgNode, E: DbgEdge> Dbg<N, E> {
    /// Convert de bruijn graph into
    /// collection of elements for cytoscape parsing
    fn to_cytoscape_elements_with_attrs_and_historys<T>(
        &self,
        node_attrs: &[NodeAttrVec],
        edge_attrs: &[EdgeAttrVec],
        node_historys: &[(String, NodeVec<DenseStorage<T>>)],
    ) -> Vec<ElementV2<N::Kmer>>
    where
        T: Into<f64> + Copy + PartialEq + Default,
    {
        let mut elements = Vec::new();

        // add history labels if node historys exists
        if node_historys.len() > 0 {
            let labels = node_historys
                .iter()
                .map(|(label, _)| label.clone())
                .collect();
            elements.push(ElementV2::HistoryLabels { labels })
        }

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
            // historys of the node
            let history = node_historys
                .iter()
                .map(|(_, vec)| vec[node].into())
                .collect();
            elements.push(ElementV2::Node {
                id: node,
                label: Some(weight.kmer().clone()),
                attrs,
                history,
            });
        }

        // add edges
        let has_copy_nums = self.is_edge_copy_nums_assigned();
        for (edge, s, t, weight) in self.edges() {
            // attrs for the edge
            let mut attrs = vec![];
            if has_copy_nums {
                attrs.push(EdgeAttr::TrueCopyNum(weight.copy_num().unwrap()));
            }
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
    pub fn to_cytoscape_with_attrs_and_historys<T>(
        &self,
        node_attrs: &[NodeAttrVec],
        edge_attrs: &[EdgeAttrVec],
        node_historys: &[(String, NodeVec<DenseStorage<T>>)],
    ) -> String
    where
        T: Into<f64> + Copy + PartialEq + Default,
    {
        let elements = self.to_cytoscape_elements_with_attrs_and_historys(
            node_attrs,
            edge_attrs,
            node_historys,
        );
        serde_json::to_string(&elements).unwrap()
    }
    //
    // wrappers of to_cytoscape_elements_with_attrs_and_historys and
    // to_cytoscape_with_attrs_and_historys.
    //
    /// Convert de bruijn graph into
    /// collection of elements for cytoscape parsing
    fn to_cytoscape_elements_with_attrs(
        &self,
        node_attrs: &[NodeAttrVec],
        edge_attrs: &[EdgeAttrVec],
    ) -> Vec<ElementV2<N::Kmer>> {
        self.to_cytoscape_elements_with_attrs_and_historys::<f64>(node_attrs, edge_attrs, &[])
    }
    /// Convert de bruijn graph into
    /// cytoscape-parsable JSON string
    /// with node/edge attributes
    pub fn to_cytoscape_with_attrs(
        &self,
        node_attrs: &[NodeAttrVec],
        edge_attrs: &[EdgeAttrVec],
    ) -> String {
        self.to_cytoscape_with_attrs_and_historys::<f64>(node_attrs, edge_attrs, &[])
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
        self.to_cytoscape_with_attrs(&[], &[])
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
    use crate::dbg::impls::SimpleDbg;
    use crate::dbg::mocks::{mock_intersection_small, mock_rep, mock_simple};
    use crate::kmer::veckmer::VecKmer;
    use crate::vector::{DenseStorage, EdgeVec, NodeVec};

    #[test]
    fn dbg_serialize() {
        let dbg = mock_simple();
        println!("{}", dbg);
        let dbg2: SimpleDbg<VecKmer> = SimpleDbg::from_str(&dbg.to_string()).unwrap();
        println!("{}", dbg2);
        assert_eq!(dbg.to_string(), dbg2.to_string());
        assert_eq!(dbg.to_string(), "4,L:AAAGCTTGATT");
        assert_eq!(dbg2.to_string(), "4,L:AAAGCTTGATT");

        let dbg = mock_rep();
        println!("{}", dbg);
        let dbg2: SimpleDbg<VecKmer> = SimpleDbg::from_str(&dbg.to_string()).unwrap();
        println!("{}", dbg2);
        assert_eq!(dbg.to_string(), "4,L:CCCCCCCCCCCCCC,L:AAAAAAAAAAAAA");
        assert_eq!(dbg2.to_string(), "4,L:CCCCCCCCCCCCCC,L:AAAAAAAAAAAAA");
    }

    #[test]
    fn dbg_mock_simple_cytoscape_with_history() {
        let mut dbg = mock_intersection_small();
        let mut v1: NodeVec<DenseStorage<_>> = NodeVec::new(dbg.n_nodes(), 0.0);
        v1[ni(1)] = 10.0;
        let mut v2: NodeVec<DenseStorage<_>> = NodeVec::new(dbg.n_nodes(), 20.0);
        v2[ni(2)] = 10.0;
        // println!("n={}", dbg.n_nodes());
        // println!("v1={}", v1);
        // println!("v2={}", v2);
        let json = dbg.to_cytoscape_with_attrs_and_historys(
            &[],
            &[],
            &[("v1".to_string(), v1), ("v2".to_string(), v2)],
        );
        println!("{}", json);
        let json_true = r#"[{"group":"history_labels","data":{"labels":["v1","v2"]}},{"group":"nodes","data":{"id":"n0","label":"ATAG","attrs":[{"type":"copy_num","value":1}],"history":[0.0,20.0]}},{"group":"nodes","data":{"id":"n1","label":"GCTn","attrs":[{"type":"copy_num","value":1}],"history":[10.0,20.0]}},{"group":"nodes","data":{"id":"n2","label":"TAAG","attrs":[{"type":"copy_num","value":1}],"history":[0.0,10.0]}},{"group":"nodes","data":{"id":"n3","label":"nnnT","attrs":[{"type":"copy_num","value":1}],"history":[0.0,20.0]}},{"group":"nodes","data":{"id":"n4","label":"nnAT","attrs":[{"type":"copy_num","value":1}],"history":[0.0,20.0]}},{"group":"nodes","data":{"id":"n5","label":"AGCC","attrs":[{"type":"copy_num","value":1}],"history":[0.0,20.0]}},{"group":"nodes","data":{"id":"n6","label":"Cnnn","attrs":[{"type":"copy_num","value":1}],"history":[0.0,20.0]}},{"group":"nodes","data":{"id":"n7","label":"GCCn","attrs":[{"type":"copy_num","value":1}],"history":[0.0,20.0]}},{"group":"nodes","data":{"id":"n8","label":"Tnnn","attrs":[{"type":"copy_num","value":1}],"history":[0.0,20.0]}},{"group":"nodes","data":{"id":"n9","label":"nnnA","attrs":[{"type":"copy_num","value":1}],"history":[0.0,20.0]}},{"group":"nodes","data":{"id":"n10","label":"nTAA","attrs":[{"type":"copy_num","value":1}],"history":[0.0,20.0]}},{"group":"nodes","data":{"id":"n11","label":"AAGC","attrs":[{"type":"copy_num","value":1}],"history":[0.0,20.0]}},{"group":"nodes","data":{"id":"n12","label":"AGCT","attrs":[{"type":"copy_num","value":1}],"history":[0.0,20.0]}},{"group":"nodes","data":{"id":"n13","label":"nATA","attrs":[{"type":"copy_num","value":1}],"history":[0.0,20.0]}},{"group":"nodes","data":{"id":"n14","label":"CTnn","attrs":[{"type":"copy_num","value":1}],"history":[0.0,20.0]}},{"group":"nodes","data":{"id":"n15","label":"nnTA","attrs":[{"type":"copy_num","value":1}],"history":[0.0,20.0]}},{"group":"nodes","data":{"id":"n16","label":"CCnn","attrs":[{"type":"copy_num","value":1}],"history":[0.0,20.0]}},{"group":"nodes","data":{"id":"n17","label":"TAGC","attrs":[{"type":"copy_num","value":1}],"history":[0.0,20.0]}},{"group":"edges","data":{"id":"e0","source":"n0","target":"n17","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e1","source":"n1","target":"n14","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e2","source":"n2","target":"n11","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e3","source":"n3","target":"n15","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e4","source":"n4","target":"n13","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e5","source":"n5","target":"n7","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e6","source":"n6","target":"n9","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e7","source":"n6","target":"n3","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e8","source":"n7","target":"n16","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e9","source":"n8","target":"n9","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e10","source":"n8","target":"n3","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e11","source":"n9","target":"n4","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e12","source":"n10","target":"n2","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e13","source":"n11","target":"n5","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e14","source":"n11","target":"n12","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e15","source":"n12","target":"n1","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e16","source":"n13","target":"n0","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e17","source":"n14","target":"n8","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e18","source":"n15","target":"n10","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e19","source":"n16","target":"n6","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e20","source":"n17","target":"n5","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e21","source":"n17","target":"n12","label":null,"attrs":[]}}]"#;
        assert_eq!(json, json_true);
    }

    #[test]
    fn dbg_mock_simple_cytoscape() {
        let dbg = mock_simple();
        // (1) with no attrs
        let json = dbg.to_cytoscape();
        println!("{}", json);
        let json_true = r#"[{"group":"nodes","data":{"id":"n0","label":"TTnn","attrs":[{"type":"copy_num","value":1}],"history":[]}},{"group":"nodes","data":{"id":"n1","label":"Tnnn","attrs":[{"type":"copy_num","value":1}],"history":[]}},{"group":"nodes","data":{"id":"n2","label":"AGCT","attrs":[{"type":"copy_num","value":1}],"history":[]}},{"group":"nodes","data":{"id":"n3","label":"GCTT","attrs":[{"type":"copy_num","value":1}],"history":[]}},{"group":"nodes","data":{"id":"n4","label":"AAAG","attrs":[{"type":"copy_num","value":1}],"history":[]}},{"group":"nodes","data":{"id":"n5","label":"CTTG","attrs":[{"type":"copy_num","value":1}],"history":[]}},{"group":"nodes","data":{"id":"n6","label":"ATTn","attrs":[{"type":"copy_num","value":1}],"history":[]}},{"group":"nodes","data":{"id":"n7","label":"TTGA","attrs":[{"type":"copy_num","value":1}],"history":[]}},{"group":"nodes","data":{"id":"n8","label":"nAAA","attrs":[{"type":"copy_num","value":1}],"history":[]}},{"group":"nodes","data":{"id":"n9","label":"nnAA","attrs":[{"type":"copy_num","value":1}],"history":[]}},{"group":"nodes","data":{"id":"n10","label":"nnnA","attrs":[{"type":"copy_num","value":1}],"history":[]}},{"group":"nodes","data":{"id":"n11","label":"AAGC","attrs":[{"type":"copy_num","value":1}],"history":[]}},{"group":"nodes","data":{"id":"n12","label":"TGAT","attrs":[{"type":"copy_num","value":1}],"history":[]}},{"group":"nodes","data":{"id":"n13","label":"GATT","attrs":[{"type":"copy_num","value":1}],"history":[]}},{"group":"edges","data":{"id":"e0","source":"n0","target":"n1","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e1","source":"n1","target":"n10","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e2","source":"n2","target":"n3","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e3","source":"n3","target":"n5","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e4","source":"n4","target":"n11","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e5","source":"n5","target":"n7","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e6","source":"n6","target":"n0","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e7","source":"n7","target":"n12","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e8","source":"n8","target":"n4","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e9","source":"n9","target":"n8","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e10","source":"n10","target":"n9","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e11","source":"n11","target":"n2","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e12","source":"n12","target":"n13","label":null,"attrs":[]}},{"group":"edges","data":{"id":"e13","source":"n13","target":"n6","label":null,"attrs":[]}}]"#;
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
        let json_true = r#"[{"group":"nodes","data":{"id":"n0","label":"TTnn","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":100},{"type":"freq","value":0.0}],"history":[]}},{"group":"nodes","data":{"id":"n1","label":"Tnnn","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":10.1}],"history":[]}},{"group":"nodes","data":{"id":"n2","label":"AGCT","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}],"history":[]}},{"group":"nodes","data":{"id":"n3","label":"GCTT","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}],"history":[]}},{"group":"nodes","data":{"id":"n4","label":"AAAG","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}],"history":[]}},{"group":"nodes","data":{"id":"n5","label":"CTTG","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}],"history":[]}},{"group":"nodes","data":{"id":"n6","label":"ATTn","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}],"history":[]}},{"group":"nodes","data":{"id":"n7","label":"TTGA","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}],"history":[]}},{"group":"nodes","data":{"id":"n8","label":"nAAA","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}],"history":[]}},{"group":"nodes","data":{"id":"n9","label":"nnAA","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}],"history":[]}},{"group":"nodes","data":{"id":"n10","label":"nnnA","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}],"history":[]}},{"group":"nodes","data":{"id":"n11","label":"AAGC","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}],"history":[]}},{"group":"nodes","data":{"id":"n12","label":"TGAT","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}],"history":[]}},{"group":"nodes","data":{"id":"n13","label":"GATT","attrs":[{"type":"copy_num","value":1},{"type":"copy_num","value":1},{"type":"freq","value":0.0}],"history":[]}},{"group":"edges","data":{"id":"e0","source":"n0","target":"n1","label":null,"attrs":[{"type":"true_copy_num","value":100}]}},{"group":"edges","data":{"id":"e1","source":"n1","target":"n10","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e2","source":"n2","target":"n3","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e3","source":"n3","target":"n5","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e4","source":"n4","target":"n11","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e5","source":"n5","target":"n7","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e6","source":"n6","target":"n0","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e7","source":"n7","target":"n12","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e8","source":"n8","target":"n4","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e9","source":"n9","target":"n8","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e10","source":"n10","target":"n9","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e11","source":"n11","target":"n2","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e12","source":"n12","target":"n13","label":null,"attrs":[{"type":"true_copy_num","value":1}]}},{"group":"edges","data":{"id":"e13","source":"n13","target":"n6","label":null,"attrs":[{"type":"true_copy_num","value":1}]}}]"#;
        assert_eq!(json, json_true);
    }
}
