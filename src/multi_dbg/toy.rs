//!
//! Toy example MultiDbg for debug
//!
use super::*;
use petgraph::graph::{DefaultIx, DiGraph, EdgeIndex, Graph, NodeIndex};

///
/// Circular `ATCG` k=4
///
/// ```text
/// k=4
///           GATC x1
///      ┌───┐      ┌───┐
///      │ATC│◄─────┤GAT│
///      └─┬─┘      └─▲─┘
///        │          │
/// ATCG x1│          │CGAT x1
///        │          │
///      ┌─▼─┐      ┌─┴─┐
///      │TCG├─────►│CGA│
///      └───┘      └───┘
///           TCGA x1
/// ```
///
pub fn circular() -> MultiDbg {
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
/// ```text
/// k=4
///         nnTA   nTAG   TAGG
///          1x     1x     1x
///     ┌───┐  ┌───┐  ┌───┐  ┌───┐
///     │nnT├─►│nTA├─►│TAG├─►│AGG│
///     └─▲─┘  └───┘  └───┘  └─┬─┘
/// nnnT  │                    │  AGGC
///  1x   │                    │   1x
///     ┌─┴─┐  ┌───┐  ┌───┐  ┌─▼─┐
///     │nnn│◄─┤Cnn│◄─┤GCn│◄─┤GGC│
///     └───┘  └───┘  └───┘  └───┘
///         Cnnn   GCnn   GGCn
///          1x     1x     1x
/// ```
pub fn linear() -> MultiDbg {
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
    let e_ggcn = full.add_edge(v_ggc, v_gcn, MultiFullEdge::new(b'n', 1));
    let e_gcnn = full.add_edge(v_gcn, v_cnn, MultiFullEdge::new(b'n', 1));
    let e_cnnn = full.add_edge(v_cnn, v_nnn, MultiFullEdge::new(b'n', 1));

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
/// ```text
/// k=4
///              ┌───┐
///   ┌──────┬───┤ATC│◄──┬──────┐
///   │      │   └───┘   │      │
/// ┌─▼─┐  ┌─▼─┐       ┌─┴─┐  ┌─┴─┐
/// │TCA│  │TCC│       │GAT│  │TAT│
/// └─┬─┘  └─┬─┘       └─▲─┘  └─▲─┘
///   │      │           │      │
/// ┌─▼─┐  ┌─▼─┐       ┌─┴─┐  ┌─┴─┐
/// │CAn│  │CCn│       │nGA│  │nTA│
/// └─┬─┘  └─┬─┘       └─▲─┘  └─▲─┘
///   │      │           │      │
/// ┌─▼─┐  ┌─▼─┐       ┌─┴─┐  ┌─┴─┐
/// │Ann│  │Cnn│       │nnG│  │nnT│
/// └─┬─┘  └─┬─┘       └─▲─┘  └─▲─┘
///   │      │   ┌───┐   │      │
///   └──────┴──►│nnn├───┴──────┘
///              └───┘
/// ```
///
pub fn intersection() -> MultiDbg {
    let mut full = DiGraph::new();
    // 2x nodes
    let v_nnn = full.add_node(MultiFullNode::new(true));
    let v_atc = full.add_node(MultiFullNode::new(false));
    // 1x nodes
    // GATCC path
    let v_nng = full.add_node(MultiFullNode::new(false));
    let v_nga = full.add_node(MultiFullNode::new(false));
    let v_gat = full.add_node(MultiFullNode::new(false));
    let v_tcc = full.add_node(MultiFullNode::new(false));
    let v_ccn = full.add_node(MultiFullNode::new(false));
    let v_cnn = full.add_node(MultiFullNode::new(false));
    // TATCA path
    let v_nnt = full.add_node(MultiFullNode::new(false));
    let v_nta = full.add_node(MultiFullNode::new(false));
    let v_tat = full.add_node(MultiFullNode::new(false));
    let v_tca = full.add_node(MultiFullNode::new(false));
    let v_can = full.add_node(MultiFullNode::new(false));
    let v_ann = full.add_node(MultiFullNode::new(false));
    // edges
    // GATCC path
    let e_nnng = full.add_edge(v_nnn, v_nng, MultiFullEdge::new(b'G', 1));
    let e_nnga = full.add_edge(v_nng, v_nga, MultiFullEdge::new(b'A', 1));
    let e_ngat = full.add_edge(v_nga, v_gat, MultiFullEdge::new(b'T', 1));
    let e_gatc = full.add_edge(v_gat, v_atc, MultiFullEdge::new(b'C', 1));
    let e_atcc = full.add_edge(v_atc, v_tcc, MultiFullEdge::new(b'C', 1));
    let e_tccn = full.add_edge(v_tcc, v_ccn, MultiFullEdge::new(b'n', 1));
    let e_ccnn = full.add_edge(v_ccn, v_cnn, MultiFullEdge::new(b'n', 1));
    let e_cnnn = full.add_edge(v_cnn, v_nnn, MultiFullEdge::new(b'n', 1));
    // TATCA path
    let e_nnnt = full.add_edge(v_nnn, v_nnt, MultiFullEdge::new(b'T', 1));
    let e_nnta = full.add_edge(v_nnt, v_nta, MultiFullEdge::new(b'A', 1));
    let e_ntat = full.add_edge(v_nta, v_tat, MultiFullEdge::new(b'T', 1));
    let e_tatc = full.add_edge(v_tat, v_atc, MultiFullEdge::new(b'C', 1));
    let e_atca = full.add_edge(v_atc, v_tca, MultiFullEdge::new(b'A', 1));
    let e_tcan = full.add_edge(v_tca, v_can, MultiFullEdge::new(b'n', 1));
    let e_cann = full.add_edge(v_can, v_ann, MultiFullEdge::new(b'n', 1));
    let e_annn = full.add_edge(v_ann, v_nnn, MultiFullEdge::new(b'n', 1));

    let compact = MultiDbg::construct_compact_from_full(&full);
    MultiDbg {
        k: 4,
        full,
        compact,
    }
}

///
/// Linear `CTAAAAAAGC` k=4
///
/// `AAAA` is 3x
///
/// ```text
/// k=4
///
///      ┌───┐
///   ┌──┤nnn│◄─┐
///   │  └───┘  │
/// ┌─▼─┐     ┌─┴─┐
/// │nnC│     │Cnn│
/// └─┬─┘     └─▲─┘
///   │         │
/// ┌─▼─┐     ┌─┴─┐
/// │nCT│     │GCn│
/// └─┬─┘     └─▲─┘
///   │         │
/// ┌─▼─┐     ┌─┴─┐
/// │CTA│     │AGC│
/// └─┬─┘     └─▲─┘
///   │         │
/// ┌─▼─┐     ┌─┴─┐
/// │TAA│     │AAG│
/// └─┬─┘     └─▲─┘
///   │  ┌───┐  │
///   └─►│   ├──┘
///      │AAA│
///   ┌─►│   ├──┐
///   │  └───┘  │
///   │         │
///   └─────────┘
/// ```
///
pub fn selfloop() -> MultiDbg {
    let mut full = DiGraph::new();
    let v_nnn = full.add_node(MultiFullNode::new(true));
    let v_nnc = full.add_node(MultiFullNode::new(false));
    let v_nct = full.add_node(MultiFullNode::new(false));
    let v_cta = full.add_node(MultiFullNode::new(false));
    let v_taa = full.add_node(MultiFullNode::new(false));
    let v_aaa = full.add_node(MultiFullNode::new(false));
    let v_aag = full.add_node(MultiFullNode::new(false));
    let v_agc = full.add_node(MultiFullNode::new(false));
    let v_gcn = full.add_node(MultiFullNode::new(false));
    let v_cnn = full.add_node(MultiFullNode::new(false));

    let e_nnnc = full.add_edge(v_nnn, v_nnc, MultiFullEdge::new(b'C', 1));
    let e_nnct = full.add_edge(v_nnc, v_nct, MultiFullEdge::new(b'T', 1));
    let e_ncta = full.add_edge(v_nct, v_cta, MultiFullEdge::new(b'A', 1));
    let e_ctaa = full.add_edge(v_cta, v_taa, MultiFullEdge::new(b'A', 1));
    let e_taaa = full.add_edge(v_taa, v_aaa, MultiFullEdge::new(b'A', 1));
    let e_aaaa = full.add_edge(v_aaa, v_aaa, MultiFullEdge::new(b'A', 3));
    let e_aaag = full.add_edge(v_aaa, v_aag, MultiFullEdge::new(b'G', 1));
    let e_aagc = full.add_edge(v_aag, v_agc, MultiFullEdge::new(b'C', 1));
    let e_agcn = full.add_edge(v_agc, v_gcn, MultiFullEdge::new(b'n', 1));
    let e_gcnn = full.add_edge(v_gcn, v_cnn, MultiFullEdge::new(b'n', 1));
    let e_cnnn = full.add_edge(v_cnn, v_nnn, MultiFullEdge::new(b'n', 1));

    let compact = MultiDbg::construct_compact_from_full(&full);
    MultiDbg {
        k: 4,
        full,
        compact,
    }
}

///
/// ```text
/// Linear k=4
///
/// `TCCCAGCAGCAGCAGGAA`
///     -->-->-->-->
/// ```
///
/// has repeat `CAG`x4
///
/// ```text
/// k=4
///       ┌───┐  ┌───┐  ┌───┐  ┌───┐  ┌───┐          ┌───┐
///   ┌───┤Ann◄──┤AAn◄──┤GAA◄──┤GGA◄──┤AGG◄───┐   ┌──►AGC├──┐
///   │   └───┘  └───┘  └───┘  └───┘  └───┘   │   │  └───┘  │
/// ┌─▼─┐                                   ┌─┴───┴─┐       │
/// │nnn│                                   │  CAG  │       │
/// └─┬─┘                                   └─▲───▲─┘       │
///   │   ┌───┐  ┌───┐  ┌───┐  ┌───┐  ┌───┐   │   │  ┌───┐  │
///   └───►nnT├──►nTC├──►TCC├──►CCC├──►CCA├───┘   └──┤GCA◄──┘
///       └───┘  └───┘  └───┘  └───┘  └───┘          └───┘
/// ```
///
pub fn repeat() -> MultiDbg {
    let mut full = DiGraph::new();
    let v_nnn = full.add_node(MultiFullNode::new(true));
    // first half
    let v_nnt = full.add_node(MultiFullNode::new(false));
    let v_ntc = full.add_node(MultiFullNode::new(false));
    let v_tcc = full.add_node(MultiFullNode::new(false));
    let v_ccc = full.add_node(MultiFullNode::new(false));
    let v_cca = full.add_node(MultiFullNode::new(false));
    // repeat unit
    let v_cag = full.add_node(MultiFullNode::new(false));
    let v_agc = full.add_node(MultiFullNode::new(false));
    let v_gca = full.add_node(MultiFullNode::new(false));
    // latter half
    let v_agg = full.add_node(MultiFullNode::new(false));
    let v_gga = full.add_node(MultiFullNode::new(false));
    let v_gaa = full.add_node(MultiFullNode::new(false));
    let v_aan = full.add_node(MultiFullNode::new(false));
    let v_ann = full.add_node(MultiFullNode::new(false));
    // edges
    let e_nnnt = full.add_edge(v_nnn, v_nnt, MultiFullEdge::new(b'T', 1));
    let e_nntc = full.add_edge(v_nnt, v_ntc, MultiFullEdge::new(b'C', 1));
    let e_ntcc = full.add_edge(v_ntc, v_tcc, MultiFullEdge::new(b'C', 1));
    let e_tccc = full.add_edge(v_tcc, v_ccc, MultiFullEdge::new(b'C', 1));
    let e_ccca = full.add_edge(v_ccc, v_cca, MultiFullEdge::new(b'A', 1));
    let e_ccag = full.add_edge(v_cca, v_cag, MultiFullEdge::new(b'G', 1));

    let e_cagc = full.add_edge(v_cag, v_agc, MultiFullEdge::new(b'C', 3));
    let e_agca = full.add_edge(v_agc, v_gca, MultiFullEdge::new(b'A', 3));
    let e_gcag = full.add_edge(v_gca, v_cag, MultiFullEdge::new(b'G', 3));

    let e_cagg = full.add_edge(v_cag, v_agg, MultiFullEdge::new(b'G', 1));
    let e_agga = full.add_edge(v_agg, v_gga, MultiFullEdge::new(b'A', 1));
    let e_ggaa = full.add_edge(v_gga, v_gaa, MultiFullEdge::new(b'A', 1));
    let e_gaan = full.add_edge(v_gaa, v_aan, MultiFullEdge::new(b'n', 1));
    let e_aann = full.add_edge(v_aan, v_ann, MultiFullEdge::new(b'n', 1));
    let e_annn = full.add_edge(v_ann, v_nnn, MultiFullEdge::new(b'n', 1));

    let compact = MultiDbg::construct_compact_from_full(&full);
    MultiDbg {
        k: 4,
        full,
        compact,
    }
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{ei, ni};

    #[test]
    fn test_circular() {
        let dbg = linear();
        dbg.show_graph_with_kmer();
        let p = dbg.get_euler_circuit();
        for s in dbg.to_styled_seqs() {
            println!("{}", s);
        }

        let dbg_ext = dbg.to_kp1_dbg();
        dbg_ext.show_graph_with_kmer();
    }
}
