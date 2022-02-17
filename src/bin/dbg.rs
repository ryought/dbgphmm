use dbgphmm::dbg::hashdbg_v2::HashDbg;
use dbgphmm::dbg::impls::SimpleDbg;
use dbgphmm::hmmv2::mocks::mock_linear_random_phmm;
use dbgphmm::hmmv2::params::PHMMParams;
use dbgphmm::kmer::veckmer::VecKmer;

fn main() {
    let phmm = mock_linear_random_phmm(100, 0, PHMMParams::default());
    println!("{}", phmm);
}
