//!
//! Toy example MultiDbg for debug
//!
use super::*;
use petgraph::graph::{DefaultIx, DiGraph, EdgeIndex, Graph, NodeIndex};

///
/// Circular `ATCG` k=4
///
fn circular() -> MultiDbg {
    let mut full = DiGraph::new();
    let v_atc = full.add_node(MultiFullNode::new(false));
    let v_tcg = full.add_node(MultiFullNode::new(false));
    let v_cga = full.add_node(MultiFullNode::new(false));
    let v_gat = full.add_node(MultiFullNode::new(false));
    let e_atcg = full.add_edge(v_atc, v_tcg, MultiFullEdge::new(b'G', 1));
    let e_tcga = full.add_edge(v_tcg, v_cga, MultiFullEdge::new(b'A', 1));
    let e_cgat = full.add_edge(v_cga, v_gat, MultiFullEdge::new(b'T', 1));
    let e_gatc = full.add_edge(v_gat, v_atc, MultiFullEdge::new(b'C', 1));

    let compact = MultiDbg::construct_compact_from_full(&full);
    MultiDbg {
        k: 4,
        full,
        compact,
    }
}

///
/// Linear `TAGGC` k=4
///
fn linear() -> MultiDbg {
    let mut full = DiGraph::new();
    let v_nnt = full.add_node(MultiFullNode::new(false));
    let v_nta = full.add_node(MultiFullNode::new(false));
    let v_tag = full.add_node(MultiFullNode::new(false));
    let v_agg = full.add_node(MultiFullNode::new(false));
    let v_ggc = full.add_node(MultiFullNode::new(false));
    let v_gcn = full.add_node(MultiFullNode::new(false));
    let v_cnn = full.add_node(MultiFullNode::new(false));
    let v_nnn = full.add_node(MultiFullNode::new(true));
    let e_nnnt = full.add_edge(v_nnn, v_nnt, MultiFullEdge::new(b'T', 1));
    let e_nnta = full.add_edge(v_nnt, v_nta, MultiFullEdge::new(b'A', 1));
    let e_ntag = full.add_edge(v_nta, v_tag, MultiFullEdge::new(b'G', 1));
    let e_tagg = full.add_edge(v_tag, v_agg, MultiFullEdge::new(b'G', 1));
    let e_aggc = full.add_edge(v_agg, v_ggc, MultiFullEdge::new(b'C', 1));
    let e_ggcn = full.add_edge(v_ggc, v_gcn, MultiFullEdge::new(b'N', 1));
    let e_gcnn = full.add_edge(v_gcn, v_cnn, MultiFullEdge::new(b'N', 1));
    let e_cnnn = full.add_edge(v_cnn, v_nnn, MultiFullEdge::new(b'N', 1));

    let compact = MultiDbg::construct_compact_from_full(&full);
    MultiDbg {
        k: 4,
        full,
        compact,
    }
}

///
/// Linear `GATCC` and `TATCA` k=4
/// Have an intersection at node `ATC`
///
fn intersection() -> MultiDbg {
    unimplemented!();
}

fn selfloop() -> MultiDbg {
    unimplemented!();
}

fn repeat() -> MultiDbg {
    unimplemented!();
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_circular() {
        // let dbg = circular();
        let dbg = linear();
        dbg.show_graph_with_kmer();
    }
}
