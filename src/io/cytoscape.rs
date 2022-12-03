use crate::common::{CopyNum, Freq};
use crate::kmer::kmer::{Kmer, KmerLike};
use crate::vector::{DenseStorage, EdgeVec, NodeVec};
use petgraph::graph::{EdgeIndex, NodeIndex};
use serde::Serialize;
use serde_with::{serde_as, DisplayFromStr};

//
// V1
//

#[serde_as]
#[derive(Debug, Clone, Serialize)]
#[serde(tag = "group", content = "data")]
pub enum Element {
    #[serde(rename = "nodes")]
    /// Node element of cytoscape
    Node {
        /// Node id
        id: usize,
        /// Node label
        #[serde_as(as = "DisplayFromStr")]
        label: Kmer,
    },
    #[serde(rename = "edges")]
    /// Edge element of cytoscape
    Edge {
        id: usize,
        source: usize,
        target: usize,
        #[serde_as(as = "DisplayFromStr")]
        label: Kmer,
        widths: Vec<u32>,
        /// TODO?
        true_width: Option<u32>,
    },
}

//
// V2
//

#[derive(Debug, Clone, Serialize)]
#[serde(tag = "type", content = "value", rename_all = "snake_case")]
pub enum NodeAttr {
    CopyNum(CopyNum),
    Freq(Freq),
}

#[derive(Debug, Clone, Serialize)]
#[serde(tag = "type", content = "value", rename_all = "snake_case")]
pub enum EdgeAttr {
    TrueCopyNum(CopyNum),
    Freq(Freq),
    CopyNums(Vec<CopyNum>),
}

#[derive(Debug, Clone)]
/// struct for storing NodeAttr for all nodes
pub enum NodeAttrVec {
    CopyNum(NodeVec<DenseStorage<CopyNum>>),
    Freq(NodeVec<DenseStorage<Freq>>),
}

#[derive(Debug, Clone)]
/// struct for storing EdgeAttr for all edges
pub enum EdgeAttrVec {
    TrueCopyNum(EdgeVec<DenseStorage<CopyNum>>),
    Freq(EdgeVec<DenseStorage<Freq>>),
}

#[serde_as]
#[derive(Debug, Clone, Serialize)]
#[serde(tag = "group", content = "data")]
pub enum ElementV2<K: KmerLike> {
    #[serde(rename = "history_labels")]
    /// Custom group name `history_labels`
    HistoryLabels {
        /// a list of labels corresponds to the history of Node
        labels: Vec<String>,
    },
    #[serde(rename = "nodes")]
    /// Node element of cytoscape
    Node {
        /// Node id
        #[serde(serialize_with = "serialize_node")]
        id: NodeIndex,
        /// Node label
        #[serde_as(as = "Option<DisplayFromStr>")]
        label: Option<K>,
        /// Node additional info
        #[serde(skip_serializing_if = "Option::is_none")]
        info: Option<String>,
        attrs: Vec<NodeAttr>,
        history: Vec<f64>,
        copy_num: CopyNum,
        /// Node expected copy num
        #[serde(skip_serializing_if = "Option::is_none")]
        copy_num_expected: Option<f64>,
    },
    #[serde(rename = "edges")]
    /// Edge element of cytoscape
    Edge {
        #[serde(serialize_with = "serialize_edge")]
        id: EdgeIndex,
        #[serde(serialize_with = "serialize_node")]
        source: NodeIndex,
        #[serde(serialize_with = "serialize_node")]
        target: NodeIndex,
        #[serde_as(as = "Option<DisplayFromStr>")]
        label: Option<K>,
        attrs: Vec<EdgeAttr>,
    },
}

fn serialize_node<S: serde::Serializer>(node: &NodeIndex, s: S) -> Result<S::Ok, S::Error> {
    s.serialize_str(&format!("n{}", node.index()))
}

fn serialize_edge<S: serde::Serializer>(edge: &EdgeIndex, s: S) -> Result<S::Ok, S::Error> {
    s.serialize_str(&format!("e{}", edge.index()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cytoscape_ser_test() {
        let mut elements = Vec::new();
        elements.push(Element::Node {
            id: 0,
            label: Kmer::from_bases(b"ATCGA"),
        });
        elements.push(Element::Edge {
            id: 0,
            source: 0,
            target: 1,
            label: Kmer::from_bases(b"ATCGAT"),
            widths: vec![0, 1, 2],
            true_width: None,
        });
        let serialized = serde_json::to_string_pretty(&elements).unwrap();
        println!("{}", serialized);
        let answer = r#"[
  {
    "group": "nodes",
    "data": {
      "id": 0,
      "label": "ATCGA"
    }
  },
  {
    "group": "edges",
    "data": {
      "id": 0,
      "source": 0,
      "target": 1,
      "label": "ATCGAT",
      "widths": [
        0,
        1,
        2
      ],
      "true_width": null
    }
  }
]"#;
        assert_eq!(serialized, answer);
    }

    #[test]
    fn cytoscape_v2_serialize_test() {
        let mut elements = Vec::new();
        elements.push(ElementV2::Node {
            id: NodeIndex::new(0),
            label: Some(Kmer::from_bases(b"ATCGA")),
            info: None,
            attrs: vec![NodeAttr::CopyNum(10)],
            history: vec![],
            copy_num: 10,
            copy_num_expected: None,
        });
        elements.push(ElementV2::Edge {
            id: EdgeIndex::new(0),
            source: NodeIndex::new(0),
            target: NodeIndex::new(1),
            label: Some(Kmer::from_bases(b"ATCGAT")),
            attrs: vec![EdgeAttr::TrueCopyNum(10), EdgeAttr::CopyNums(vec![10, 20])],
        });
        let serialized = serde_json::to_string_pretty(&elements).unwrap();
        println!("{}", serialized);
        let answer = r#"[
  {
    "group": "nodes",
    "data": {
      "id": "n0",
      "label": "ATCGA",
      "attrs": [
        {
          "type": "copy_num",
          "value": 10
        }
      ],
      "history": [],
      "copy_num": 10
    }
  },
  {
    "group": "edges",
    "data": {
      "id": "e0",
      "source": "n0",
      "target": "n1",
      "label": "ATCGAT",
      "attrs": [
        {
          "type": "true_copy_num",
          "value": 10
        },
        {
          "type": "copy_nums",
          "value": [
            10,
            20
          ]
        }
      ]
    }
  }
]"#;
        assert_eq!(serialized, answer);
    }
}
