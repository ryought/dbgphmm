use dbgphmm::dbg::hashdbg_v2::HashDbg;
use dbgphmm::dbg::impls::SimpleDbg;
use dbgphmm::graph::mocks::mock_crossing;
use dbgphmm::graph::seq_graph::SeqGraph;
use dbgphmm::hmmv2::mocks::mock_linear_random_phmm;
use dbgphmm::hmmv2::params::PHMMParams;
use dbgphmm::kmer::veckmer::VecKmer;

fn main() {
    /*
    let gg = mock_crossing(true);
    println!("{}", gg);
    let phmm = mock_linear_random_phmm(100, 0, PHMMParams::default());
    */
    let g1 = mock_crossing(false)
        .to_seq_graph()
        .to_phmm(PHMMParams::default());
    println!("{}", g1);
}
