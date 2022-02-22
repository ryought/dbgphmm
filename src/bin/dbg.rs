use dbgphmm::dbg::hashdbg_v2::HashDbg;
use dbgphmm::dbg::impls::SimpleDbg;
use dbgphmm::dbg::mocks::{mock_random, mock_simple};
use dbgphmm::graph::mocks::mock_crossing;
use dbgphmm::graph::seq_graph::SeqGraph;
use dbgphmm::hmmv2::mocks::mock_linear_random_phmm;
use dbgphmm::hmmv2::params::PHMMParams;
use dbgphmm::io::cytoscape::NodeAttrVec;
use dbgphmm::kmer::veckmer::VecKmer;

fn main() {
    /*
    let gg = mock_crossing(true);
    println!("{}", gg);
    let phmm = mock_linear_random_phmm(100, 0, PHMMParams::default());
    let g1 = mock_crossing(false)
        .to_seq_graph()
        .to_phmm(PHMMParams::default());
    println!("{}", g1);
    */
    // let dbg = mock_random(8, 1000);
    let dbg = mock_simple();
    // let json = dbg.to_edbg().to_cytoscape();
    // println!("{}", json);
    let phmm = dbg.to_phmm(PHMMParams::default());
    let r = b"GCTTGA";
    let o = phmm.run(r);
    let node_freqs = o.to_node_freqs();
    let json = dbg.to_cytoscape_with_attrs(&[NodeAttrVec::Freq(node_freqs)], &[]);
    println!("{}", json);
}
