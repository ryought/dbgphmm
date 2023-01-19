use clap::Clap;
use dbgphmm::dbg::greedy::get_max_posterior_instance;
use dbgphmm::dbg::hashdbg_v2::HashDbg;
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

    for seed in 0..3 {
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

        let mut k = opts.k_init;

        // inference
        while k <= opts.k_final {
            let mut dbg: SimpleDbg<VecKmer> = SimpleDbg::from_styled_seqs(k, &genome);
            let copy_nums_true = dbg.to_node_copy_nums();
            println!("# k={} s={} n_nodes={}", k, seed, dbg.n_nodes());
            println!("# k={} s={} n_edges={}", k, seed, dbg.n_edges());
            println!(
                "# k={} s={} copy_num_stats={:?}",
                k,
                seed,
                dbg.copy_num_stats()
            );
            println!("# k={} s={} degree_stats={:?}", k, seed, dbg.degree_stats());

            let distribution =
                dbg.search_posterior_raw(&genome, param, 10, 3, genome_size, 100, |_instance| {});

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
                    instance.info_string(),
                    kmers_to_string_pretty(&missings),
                    kmers_to_string_pretty(&errors),
                );
            }

            // dbg.set_node_copy_nums(&copy_nums_true);
            dbg.set_node_copy_nums(get_max_posterior_instance(&distribution).copy_nums());
            let neighbors: Vec<_> = distribution
                .iter()
                .map(|(p, instance, _score)| (instance.copy_nums().clone(), *p))
                .collect();
            let read_count = HashDbg::from_styled_seqs(k, &genome);
            dbg.inspect_kmer_variance_with_comment(
                &neighbors,
                &copy_nums_true,
                &read_count,
                |_| format!("{}", seed),
            );

            // upgrade
            dbg = dbg.to_k_max_dbg_naive(opts.k_final);
            if k == dbg.k() {
                break;
            }
            k = dbg.k();
        }
    }
}
