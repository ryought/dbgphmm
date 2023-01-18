use clap::Clap;
use dbgphmm::common::{ni, Reads};
use dbgphmm::dbg::float::{q_score_diff_exact, CopyDensity, FloatDbg, FloatDbgEdge, FloatDbgNode};
use dbgphmm::dbg::mocks::mock_intersection_small;
use dbgphmm::dbg::{Dbg, SimpleDbg};
use dbgphmm::genome;
use dbgphmm::kmer::VecKmer;
use dbgphmm::prelude::*;
use dbgphmm::vector::{DenseStorage, NodeVec};

///
/// Generate genome de Bruijn graph for dbgviz visualization
///
#[derive(Clap, Debug)]
struct Opts {
    ///
    #[clap(short = 'k')]
    k: usize,
    ///
    #[clap(short = 'U')]
    unit_size: usize,
    ///
    #[clap(short = 'N')]
    n_unit: usize,
    ///
    #[clap(short = 'H', default_value = "0.01")]
    hap_divergence: f64,
    ///
    #[clap(short = 'P', default_value = "1")]
    n_haplotypes: usize,
    ///
    #[clap(short = 's')]
    seed: u64,
}

fn main() {
    let opts: Opts = Opts::parse();
    // data generation
    let (genome, genome_size) = genome::tandem_repeat_polyploid_with_unique_homo_ends(
        opts.unit_size,
        opts.n_unit,
        opts.seed,
        50,
        opts.n_haplotypes,
        opts.hap_divergence,
        opts.seed,
    );
    let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_styled_seqs(opts.k, &genome);
    let json = dbg.to_cytoscape();
    println!("{}", json);
}

fn main_old() {
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
