use dbgphmm::common::{ni, Reads};
use dbgphmm::dbg::float::{q_score_diff_exact, CopyDensity, FloatDbg, FloatDbgEdge, FloatDbgNode};
use dbgphmm::dbg::mocks::mock_intersection_small;
use dbgphmm::em::float::{em, em_result_to_node_historys};
use dbgphmm::prelude::*;
use dbgphmm::vector::{DenseStorage, NodeVec};

fn main() {
    let dbg = mock_intersection_small();
    let genome_size = dbg.genome_size() as CopyDensity;
    let reads = Reads::from(vec![b"ATAGCT".to_vec()]);
    let fdbg = FloatDbg::from_dbg(&dbg);
    let params = PHMMParams::zero_error();

    let result = em(&fdbg, &reads, &params, genome_size, 0.1, 10, 10);
    let historys = em_result_to_node_historys(&result);

    let json = dbg.to_cytoscape_with_attrs_and_historys(&[], &[], &historys);
    // let (edge_freqs, init_freqs, p) = e_step(&fdbg, &reads, &params);
    println!("{}", json);
}
