use clap::Parser;
use dbgphmm::{
    e2e::{generate_dataset, ReadType},
    genome,
    hmmv2::params::PHMMParams,
    multi_dbg::MultiDbg,
    utils::{check_memory_usage, timer},
};

#[derive(Parser, Debug)]
struct Opts {
    #[clap(short)]
    k: usize,
    #[clap(short = 'C')]
    coverage: usize,
    #[clap(short = 'L')]
    read_length: usize,
    #[clap(short = 'p')]
    p_error: f64,
    // genome related
    #[clap(short = 'U')]
    unit_size: usize,
    #[clap(short = 'N')]
    n_unit: usize,
    #[clap(short = 'E')]
    end_length: usize,
    #[clap(short = 'H')]
    hap_divergence: f64,
    #[clap(short = 'P')]
    n_haplotypes: usize,
    // output
    #[clap(long)]
    output_prefix: std::path::PathBuf,
    #[clap(long)]
    genome_fasta: Option<std::path::PathBuf>,
}

// -U 10000 -N 10 -E 200 -P 2 -H 0.01
// -C 20 -p 0.01 -L 20000

fn main() {
    let opts: Opts = Opts::parse();

    println!("# started_at={}", chrono::Local::now());
    println!("# opts={:?}", opts);
    check_memory_usage();

    let (genome, genome_size) = genome::tandem_repeat_polyploid_with_unique_homo_ends(
        opts.unit_size,
        opts.n_unit,
        0,
        opts.end_length,
        opts.n_haplotypes,
        opts.hap_divergence,
        0,
    );
    let param = PHMMParams::uniform(opts.p_error);
    let (dataset, t) = timer(|| {
        generate_dataset(
            genome.clone(),
            genome_size,
            0,
            opts.coverage,
            opts.read_length,
            ReadType::FragmentWithRevComp,
            param,
        )
    });
    println!("# dataset created in {}ms", t);

    // dataset.show_reads_with_genome();
    let (_, t) = timer(|| dataset.to_json_file(opts.output_prefix.with_extension("json")));
    println!("# dataset dumped in {}ms", t);

    let (mut mdbg, t) = timer(|| MultiDbg::create_draft_from_dataset(opts.k, &dataset));
    println!("# draft dbg created in {}ms", t);
    mdbg.to_gfa_file(opts.output_prefix.with_extension("gfa"));
    mdbg.to_dbg_file(opts.output_prefix.with_extension("dbg"));

    match mdbg.compact_paths_from_styled_seqs(dataset.genome()) {
        Ok(paths_true) => {
            mdbg.to_paths_file(opts.output_prefix.with_extension("paths"), &paths_true);
        }
        Err(notfound) => {
            println!("# kmer notfound {}", notfound);
        }
    };

    println!("# finished_at={}", chrono::Local::now());
}
