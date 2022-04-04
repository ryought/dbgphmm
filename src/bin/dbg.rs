use dbgphmm::common::Reads;
use dbgphmm::dbg::hashdbg_v2::HashDbg;
use dbgphmm::dbg::impls::SimpleDbg;
use dbgphmm::dbg::mocks::{mock_intersection, mock_random, mock_simple};
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

    /*
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
    */

    /*
    let mut dbg = mock_intersection();
    let (ncn, ecn) = dbg.to_copy_nums_of_seq(b"AACTAGCTT").unwrap();
    dbg.set_node_copy_nums(&ncn);
    dbg.set_edge_copy_nums(Some(&ecn));
    println!("{}", dbg.to_cytoscape());
    */

    /*
    let mut dbg = mock_intersection();
    // println!("{}", dbg);
    let (ncn, ecn) = dbg.to_copy_nums_of_seq(b"AACTAGCTT").unwrap();
    dbg.set_node_copy_nums(&ncn);
    dbg.set_edge_copy_nums(Some(&ecn));
    // println!("{}", dbg);
    let dbg2 = dbg.to_kp1_dbg();
    println!("{}", dbg2.to_cytoscape());
    // println!("{}", dbg2);
    */

    let reads = Reads {
        reads: vec![b"AACTAGCTT".to_vec(), b"AACTAGCTT".to_vec()],
    };

    let mut dbg = mock_intersection();

    // (a)
    let phmm = dbg.to_phmm(PHMMParams::default());
    // println!("{}", phmm);
    let nf = phmm.to_node_freqs(&reads);
    // println!("{}", nf);
    /*
    println!(
        "{}",
        dbg.to_cytoscape_with_attrs(&[NodeAttrVec::Freq(nf)], &[])
    );
    */

    // changing edge copy nums
    let (ncn, ecn) = dbg.to_copy_nums_of_seq(b"CCGTAGGGC").unwrap();
    dbg.set_node_copy_nums(&ncn);
    dbg.set_edge_copy_nums(Some(&ecn));

    // (b)
    let phmm = dbg.to_phmm(PHMMParams::default());
    // println!("{}", phmm);
    let nf = phmm.to_node_freqs(&reads);
    println!("{}", nf);
    /*
    println!(
        "{}",
        dbg.to_cytoscape_with_attrs(&[NodeAttrVec::Freq(nf)], &[])
    );
    */
}
