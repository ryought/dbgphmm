//!
//! Node mapping
//!
//! If G=(V,E) is converted into G'=(V',E'), we want to convert a list of nodes W ⊆ V into W' ⊆ V'.
//!
//! Mapping from G node into G' node(s) `m(v)={v', ...}`
//!
use petgraph::graph::NodeIndex;

struct NodeMap(Vec<Vec<NodeIndex>>);

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn a() {}
}
