use clap::Clap;
use dbgphmm::dbg::{Dbg, SimpleDbg};
use dbgphmm::genome;
use dbgphmm::hmmv2::params::PHMMParams;
use dbgphmm::kmer::common::kmers_to_string_pretty;
use dbgphmm::kmer::VecKmer;

///
/// Sample posterior distribution batch tester to search for examples P(G_true) < P(G')
///
/// * Use genome full length read
///
/// * Test different k (k0 <= k < K) and multiple seeds at once
///
#[derive(Clap, Debug)]
struct Opts {
    ///
    #[clap(short = 'k', default_value = "12")]
    k_init: usize,
    ///
    #[clap(short = 'K', default_value = "40")]
    k_final: usize,
    /// Read error rate
    #[clap(short = 'p', default_value = "0.01")]
    p_error: f64,
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
}

fn main() {
    let opts: Opts = Opts::parse();
    let param = PHMMParams::uniform(opts.p_error);
    println!("# opts={:?}", opts);

    for seed in 0..5 {
        // data generation
        let (genome, genome_size) = genome::tandem_repeat_polyploid_with_unique_homo_ends(
            opts.unit_size,
            opts.n_unit,
            seed,
            50,
            opts.n_haplotypes,
            opts.hap_divergence,
            seed,
        );
        println!("# seed={}", seed);
        println!("# genome_size={}", genome_size);
        for i in 0..genome.len() {
            println!("# genome[{}]={}", i, genome[i]);
        }

        // inference
        for k in opts.k_init..=opts.k_final {
            let mut dbg: SimpleDbg<VecKmer> = SimpleDbg::from_styled_seqs(k, &genome);
            let copy_nums_true = dbg.to_node_copy_nums();
            println!("# k={} n_nodes={}", k, dbg.n_nodes());
            println!("# k={} n_edges={}", k, dbg.n_edges());
            println!("# k={} copy_num_stats={:?}", k, dbg.copy_num_stats());
            println!("# k={} degree_stats={:?}", k, dbg.degree_stats());

            // TODO not use dataset but use seq list
            let distribution =
                dbg.search_posterior_raw(&genome, param, 10, 3, genome_size, 100, |instance| {});

            for (p, instance, score) in distribution.iter() {
                dbg.set_node_copy_nums(instance.copy_nums());
                let ((n_missing, _n_missing_null), (n_error, _n_error_null)) =
                    dbg.inspect_kmers(&genome);
                let (missings, errors) = dbg.missing_error_kmers(&genome);
                println!(
                    "N\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    k,
                    seed,
                    p,
                    score.p_rg().to_log_value(),
                    score.p_g().to_log_value(),
                    dbg.genome_size(),
                    dbg.n_starting_kmers(),
                    instance.move_count(),
                    instance.copy_nums().dist(&copy_nums_true),
                    instance.copy_nums().max_abs_diff(&copy_nums_true),
                    n_missing,
                    n_error,
                    instance.info(),
                    kmers_to_string_pretty(&missings),
                    kmers_to_string_pretty(&errors),
                );
            }
        }
    }
}
