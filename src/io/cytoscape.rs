use crate::kmer::kmer::{Kmer, KmerLike};
use serde::Serialize;
use serde_with::{serde_as, DisplayFromStr};

#[serde_as]
#[derive(Debug, Clone, Serialize)]
#[serde(tag = "group", content = "data")]
pub enum Element {
    #[serde(rename = "nodes")]
    Node {
        id: usize,
        #[serde_as(as = "DisplayFromStr")]
        label: Kmer,
    },
    #[serde(rename = "edges")]
    Edge {
        id: usize,
        source: usize,
        target: usize,
        #[serde_as(as = "DisplayFromStr")]
        label: Kmer,
        widths: Vec<u32>,
        true_width: Option<u32>,
    },
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
    }
}
