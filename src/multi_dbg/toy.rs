#![allow(unused)]
//!
//! Toy example MultiDbg for debug
//!
//! * circular
//! * linear
//! * intersection
//! * selfloop
//! * repeat
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
    // e0-e5
    let e_nnnt = full.add_edge(v_nnn, v_nnt, MultiFullEdge::new(b'T', 1));
    let e_nntc = full.add_edge(v_nnt, v_ntc, MultiFullEdge::new(b'C', 1));
    let e_ntcc = full.add_edge(v_ntc, v_tcc, MultiFullEdge::new(b'C', 1));
    let e_tccc = full.add_edge(v_tcc, v_ccc, MultiFullEdge::new(b'C', 1));
    let e_ccca = full.add_edge(v_ccc, v_cca, MultiFullEdge::new(b'A', 1));
    let e_ccag = full.add_edge(v_cca, v_cag, MultiFullEdge::new(b'G', 1));
    // e6-e8
    let e_cagc = full.add_edge(v_cag, v_agc, MultiFullEdge::new(b'C', 3));
    let e_agca = full.add_edge(v_agc, v_gca, MultiFullEdge::new(b'A', 3));
    let e_gcag = full.add_edge(v_gca, v_cag, MultiFullEdge::new(b'G', 3));
    // e9-e14
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

///
/// ```text
/// Linear k=5
///
/// `TCCCAGCAGCAGCAGGAA`
///     -->-->-->-->
/// ```
///
/// ```text
/// k=4
///
///   1x  ┌───┐  ┌───┐  ┌───┐  ┌───┐  ┌───┐  1x   3x ┌───┐
///   ┌───┤Ann◄──┤AAn◄──┤GAA◄──┤GGA◄──┤AGG◄───┐   ┌──►AGC├──┐
///   │   └───┘  └───┘  └───┘  └───┘  └───┘   │   │  └───┘  │
/// ┌─▼─┐                                   ┌─┴───┴─┐       │
/// │nnn│                                   │  CAG  │       │ 3x
/// └─┬─┘                                   └─▲───▲─┘       │
///   │   ┌───┐  ┌───┐  ┌───┐  ┌───┐  ┌───┐   │   │  ┌───┐  │
///   └───►nnT├──►nTC├──►TCC├──►CCC├──►CCA├───┘   └──┤GCA◄──┘
///   1x  └───┘  └───┘  └───┘  └───┘  └───┘  1x   3x └───┘
///
///
///
/// k=5
///
///   1x        1x      1x      1x      1x      1x      1x        3x
/// Annnn     AAnnn   GAAnn   GGAAn   AGGAA   CAGGA   GCAGG     AGCAG
///       ┌────┐  ┌────┐  ┌────┐  ┌────┐  ┌────┐  ┌────┐  ┌────┐
///    ┌──┤Annn◄──┤AAnn◄──┤GAAn◄──┤GGAA◄──┤AGGA◄──┤CAGG◄──┤GCAG◄─────┐
///    │  └────┘  └────┘  └────┘  └────┘  └────┘  └▲───┘  └───┬┘     │
/// ┌──▼─┐                                         │          │    ┌─┴──┐
/// │nnnn│                                    CCAGG│     GCAGC│    │AGCA│
/// └──┬─┘                                      0x │       2x │    └─▲──┘
///    │  ┌────┐  ┌────┐  ┌────┐  ┌────┐  ┌────┐  ┌┴───┐  ┌───▼┐     │
///    └──►nnnT├──►nnTC├──►nTCC├──►TCCC├──►CCCA├──►CCAG├──►CAGC├─────┘
///       └────┘  └────┘  └────┘  └────┘  └────┘  └────┘  └────┘
/// nnnnT     nnnTC   nnTCC   nTCCC   TCCCA   CCCAG   CCAGC     CAGCA
///   1x        1x      1x      1x      1x      1x      1x        3x
///
///           e4                                  ┌────┐e0┌────┐
///    ┌──────────────────────────────────────────┤CAGG◄──┤GCAG◄─────┐
///    │                                          └▲───┘  └───┬┘     │
/// ┌──▼─┐                                         │          │      │
/// │nnnn│                                       e2│        e5│    e3│
/// └──┬─┘                                         │          │      │
///    │                                          ┌┴───┐  ┌───▼┐     │
///    └──────────────────────────────────────────►CCAG├──►CAGC├─────┘
///           e6                                  └────┘e1└────┘
/// ```
///
pub fn repeat_kp1() -> MultiDbg {
    let mut full = DiGraph::new();
    // nodes
    // first half
    let v_nnnt = full.add_node(MultiFullNode::new(false));
    let v_nntc = full.add_node(MultiFullNode::new(false));
    let v_ntcc = full.add_node(MultiFullNode::new(false));
    let v_tccc = full.add_node(MultiFullNode::new(false));
    let v_ccca = full.add_node(MultiFullNode::new(false));
    //
    let v_ccag = full.add_node(MultiFullNode::new(false));
    let v_cagc = full.add_node(MultiFullNode::new(false));
    let v_agca = full.add_node(MultiFullNode::new(false));
    let v_gcag = full.add_node(MultiFullNode::new(false));
    let v_cagg = full.add_node(MultiFullNode::new(false));
    // last half
    let v_agga = full.add_node(MultiFullNode::new(false));
    let v_ggaa = full.add_node(MultiFullNode::new(false));
    let v_gaan = full.add_node(MultiFullNode::new(false));
    let v_aann = full.add_node(MultiFullNode::new(false));
    let v_annn = full.add_node(MultiFullNode::new(false));
    // terminal
    let v_nnnn = full.add_node(MultiFullNode::new(true));

    // edges
    // e0-e1
    let e_annnn = full.add_edge(v_annn, v_nnnn, MultiFullEdge::new(b'n', 1));
    let e_nnnnt = full.add_edge(v_nnnn, v_nnnt, MultiFullEdge::new(b'T', 1));
    // e2-e6
    let e_nnntc = full.add_edge(v_nnnt, v_nntc, MultiFullEdge::new(b'C', 1));
    let e_nntcc = full.add_edge(v_nntc, v_ntcc, MultiFullEdge::new(b'C', 1));
    let e_ntccc = full.add_edge(v_ntcc, v_tccc, MultiFullEdge::new(b'C', 1));
    let e_tccca = full.add_edge(v_tccc, v_ccca, MultiFullEdge::new(b'A', 1));
    let e_cccag = full.add_edge(v_ccca, v_ccag, MultiFullEdge::new(b'G', 1));
    // e7-e10
    let e_gcagg = full.add_edge(v_gcag, v_cagg, MultiFullEdge::new(b'G', 1));
    let e_gcagc = full.add_edge(v_gcag, v_cagc, MultiFullEdge::new(b'C', 2));
    let e_ccagg = full.add_edge(v_ccag, v_cagg, MultiFullEdge::new(b'G', 0));
    let e_ccagc = full.add_edge(v_ccag, v_cagc, MultiFullEdge::new(b'C', 1));
    // e11-e12
    let e_cagca = full.add_edge(v_cagc, v_agca, MultiFullEdge::new(b'A', 3));
    let e_agcag = full.add_edge(v_agca, v_gcag, MultiFullEdge::new(b'G', 3));
    // e13-e17
    let e_cagga = full.add_edge(v_cagg, v_agga, MultiFullEdge::new(b'A', 1));
    let e_aggaa = full.add_edge(v_agga, v_ggaa, MultiFullEdge::new(b'A', 1));
    let e_ggaan = full.add_edge(v_ggaa, v_gaan, MultiFullEdge::new(b'n', 1));
    let e_gaann = full.add_edge(v_gaan, v_aann, MultiFullEdge::new(b'n', 1));
    let e_aannn = full.add_edge(v_aann, v_annn, MultiFullEdge::new(b'n', 1));

    let compact = MultiDbg::construct_compact_from_full(&full);
    MultiDbg {
        k: 5,
        full,
        compact,
    }
}

/// k=3 DBG
///
/// ```text
/// GGACGT GGACGT GGACGT GGACGT GGAAGT
///                                x
/// ```
///
/// |              | sequence | abundance |
/// | ------------ | -------- | --------- |
/// | unit         | `GGACGT` | 4         |
/// | mutated unit | `GGAAGT` | 1         |
///
/// ```text
///      nnG 1x       ┌──┐    Tnn 1x
///     ┌─────────────┤nn◀─────────────┐
///     │             └──┘             │
///   ┌─▼┐                            ┌┴─┐
///   │nG│                            │Tn│
///   └─┬┘                            └▲─┘
///  nGG│                              │GTn
///   1x│                              │1x
///   ┌─▼┐ TGG 4x     ┌──┐ GTG 4x     ┌┴─┐
///   │GG◀────────────┤TG◀────────────┤GT│
///   └─┬┘            └──┘            └▲─▲
///  GGA│                              │ │
///   5x│                              │ │
///   ┌─▼┐ GAC  ┌──┐ ACG  ┌──┐ CGT     │ │
///   │GA├──────▶AC├──────▶CG├─────────┘ │
///   └─┬┘  4x  └──┘  4x  └──┘  4x       │
///  GAA│                                │
///   1x│                                │
///     │  GAA  ┌──┐ AAG  ┌──┐ AGT       │
///     └───────▶AA├──────▶AG├───────────┘
///         1x  └──┘  1x  └──┘  1x
/// ```
///
pub fn one_in_n_repeat() -> MultiDbg {
    let mut full = DiGraph::new();
    // Nodes
    let v_nn = full.add_node(MultiFullNode::new(true));
    // in
    let v_ng = full.add_node(MultiFullNode::new(false));
    // major unit
    let v_gg = full.add_node(MultiFullNode::new(false));
    let v_ga = full.add_node(MultiFullNode::new(false));
    let v_ac = full.add_node(MultiFullNode::new(false));
    let v_cg = full.add_node(MultiFullNode::new(false));
    let v_gt = full.add_node(MultiFullNode::new(false));
    let v_tg = full.add_node(MultiFullNode::new(false));
    // mutated unit
    let v_aa = full.add_node(MultiFullNode::new(false));
    let v_ag = full.add_node(MultiFullNode::new(false));
    // out
    let v_tn = full.add_node(MultiFullNode::new(false));

    // Edges
    let e_nng = full.add_edge(v_nn, v_ng, MultiFullEdge::new(b'C', 1));
    let e_ngg = full.add_edge(v_ng, v_gg, MultiFullEdge::new(b'C', 1));
    // major
    let e_gga = full.add_edge(v_gg, v_ga, MultiFullEdge::new(b'A', 5));
    let e_gac = full.add_edge(v_ga, v_ac, MultiFullEdge::new(b'C', 4));
    let e_acg = full.add_edge(v_ac, v_cg, MultiFullEdge::new(b'G', 4));
    let e_cgt = full.add_edge(v_cg, v_gt, MultiFullEdge::new(b'T', 4));
    let e_gtg = full.add_edge(v_gt, v_tg, MultiFullEdge::new(b'G', 4));
    let e_tgg = full.add_edge(v_tg, v_gg, MultiFullEdge::new(b'G', 4));
    // mutated
    let e_gaa = full.add_edge(v_ga, v_aa, MultiFullEdge::new(b'A', 1));
    let e_aag = full.add_edge(v_aa, v_ag, MultiFullEdge::new(b'G', 1));
    let e_agt = full.add_edge(v_ag, v_gt, MultiFullEdge::new(b'T', 1));
    // out
    let e_gtn = full.add_edge(v_gt, v_tn, MultiFullEdge::new(b'n', 1));
    let e_tnn = full.add_edge(v_tn, v_nn, MultiFullEdge::new(b'n', 1));

    let compact = MultiDbg::construct_compact_from_full(&full);
    MultiDbg {
        k: 3,
        full,
        compact,
    }
}

/// k=3 DBG
///
/// ```text
/// Linear
///     ACTTG
///     AATTG
/// Circular
///     CAGG
/// ```
///
/// ```text
///         Gnn 2x   ┌──┐ TGn 2x
///        ┌─────────┤Gn◀────────┐
///        │         └──┘        │           ┌──┐ GCA ┌──┐
///       ┌▼─┐                 ┌─┴┐          │CA◀─────┤GC│
///       │nn│                 │TG│          └┬─┘  2x └─▲┘
///       └┬─┘                 └─▲┘        CAG│2x       │GGC 2x
///  nnA 2x│                     │TTG 2x     ┌▼─┐ AGG ┌─┴┐
///       ┌▼─┐nAC┌──┐ACT┌──┐CTT┌─┴┐          │AG├─────▶GG│
///       │nA├───▶AC├───▶CT├───▶TT│          └──┘  2x └──┘
///       └┬─┘ 1x└──┘ 1x└──┘ 1x└─▲┘
///        │                     │
///  nAA 1x│  nAA┌──┐AAT┌──┐ATT  │
///        └─────▶AA├───▶AT├─────┘
///            1x└──┘ 1x└──┘ 1x
/// ```
///
pub fn two_components() -> MultiDbg {
    let mut full = DiGraph::new();

    // Nodes
    // component1
    // common
    let v_nn = full.add_node(MultiFullNode::new(true));
    let v_na = full.add_node(MultiFullNode::new(false));
    // branch1 or branch2
    let v_tt = full.add_node(MultiFullNode::new(false));
    let v_tg = full.add_node(MultiFullNode::new(false));
    let v_gn = full.add_node(MultiFullNode::new(false));
    // branch1
    let v_ac = full.add_node(MultiFullNode::new(false));
    let v_ct = full.add_node(MultiFullNode::new(false));
    // branch2
    let v_aa = full.add_node(MultiFullNode::new(false));
    let v_at = full.add_node(MultiFullNode::new(false));
    // component2
    let v_ca = full.add_node(MultiFullNode::new(false));
    let v_ag = full.add_node(MultiFullNode::new(false));
    let v_gg = full.add_node(MultiFullNode::new(false));
    let v_gc = full.add_node(MultiFullNode::new(false));

    // Edges
    // component1
    // common
    let e_nna = full.add_edge(v_nn, v_na, MultiFullEdge::new(b'A', 2));
    // branch1
    let e_nac = full.add_edge(v_na, v_ac, MultiFullEdge::new(b'C', 1));
    let e_act = full.add_edge(v_ac, v_ct, MultiFullEdge::new(b'T', 1));
    let e_ctt = full.add_edge(v_ct, v_tt, MultiFullEdge::new(b'T', 1));
    // branch2
    let e_naa = full.add_edge(v_na, v_aa, MultiFullEdge::new(b'A', 1));
    let e_aat = full.add_edge(v_aa, v_at, MultiFullEdge::new(b'T', 1));
    let e_att = full.add_edge(v_at, v_tt, MultiFullEdge::new(b'T', 1));
    // common
    let e_ttg = full.add_edge(v_tt, v_tg, MultiFullEdge::new(b'G', 2));
    let e_tgn = full.add_edge(v_tg, v_gn, MultiFullEdge::new(b'n', 2));
    let e_gnn = full.add_edge(v_gn, v_nn, MultiFullEdge::new(b'n', 2));
    // component2
    let e_cag = full.add_edge(v_ca, v_ag, MultiFullEdge::new(b'G', 2));
    let e_agg = full.add_edge(v_ag, v_gg, MultiFullEdge::new(b'G', 2));
    let e_ggc = full.add_edge(v_gg, v_gc, MultiFullEdge::new(b'C', 2));
    let e_gca = full.add_edge(v_gc, v_ca, MultiFullEdge::new(b'A', 2));

    let compact = MultiDbg::construct_compact_from_full(&full);
    MultiDbg {
        k: 3,
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
    use crate::common::{ei, ni, StyledSequence};

    #[test]
    fn test_circular() {
        let dbg = circular();
        dbg.show_graph_with_kmer();
        let p = dbg.get_euler_circuits();

        let dbg_ext = dbg.to_kp1_dbg();
        dbg_ext.show_graph_with_kmer();
    }

    #[test]
    fn test_intersection() {
        let dbg = intersection();
        dbg.show_graph_with_kmer();

        let p = dbg.get_euler_circuits();
        assert_eq!(
            p,
            vec![
                vec![ei(0), ei(1), ei(2), ei(3), ei(12), ei(13), ei(14), ei(15)],
                vec![ei(8), ei(9), ei(10), ei(11), ei(4), ei(5), ei(6), ei(7)]
            ]
        );
        for s in dbg.to_styled_seqs() {
            println!("{}", s);
        }
        assert_eq!(
            dbg.to_styled_seqs(),
            vec![
                StyledSequence::linear(b"GATCA".to_vec()),
                StyledSequence::linear(b"TATCC".to_vec()),
            ]
        );
    }

    #[test]
    fn test_selfloop() {
        let dbg = selfloop();
        dbg.show_graph_with_kmer();

        let p = dbg.get_euler_circuits();
        assert_eq!(
            p,
            vec![vec![
                ei(0),
                ei(1),
                ei(2),
                ei(3),
                ei(4),
                ei(5),
                ei(5),
                ei(5),
                ei(6),
                ei(7),
                ei(8),
                ei(9),
                ei(10)
            ]]
        );
        for s in dbg.to_styled_seqs() {
            println!("{}", s);
        }
        assert_eq!(
            dbg.to_styled_seqs(),
            vec![StyledSequence::linear(b"CTAAAAAAGC".to_vec())]
        );
    }
}
