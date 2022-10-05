use dbgphmm::common::{ni, Reads};
use dbgphmm::dbg::float::{q_score_diff_exact, CopyDensity, FloatDbg, FloatDbgEdge, FloatDbgNode};
use dbgphmm::dbg::mocks::mock_intersection_small;
use dbgphmm::prelude::*;
use dbgphmm::vector::{DenseStorage, NodeVec};

fn main() {
    let dbg = mock_intersection_small();
    let genome_size = dbg.genome_size() as CopyDensity;
    let reads = Reads::from(vec![b"ATAGCT".to_vec()]);
    let fdbg = FloatDbg::from_dbg(&dbg);
    let params = PHMMParams::zero_error();

    let phmm = dbg.to_phmm(params.clone());
    let (nf, p) = phmm.to_node_freqs_parallel(&reads);

    // println!("{}", nf);

    let json = dbg.to_cytoscape_with_attrs_and_historys(
        &[],
        &[],
        &[
            ("node_freqs".to_string(), nf.clone()),
            ("node_freqs 2".to_string(), nf),
        ],
    );
    // let (edge_freqs, init_freqs, p) = e_step(&fdbg, &reads, &params);
    println!("{}", json);
}
