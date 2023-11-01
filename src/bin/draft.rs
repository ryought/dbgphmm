use clap::Parser;
use dbgphmm::{
    e2e::{generate_dataset, ReadType},
    genome,
    hashdbg::HashDbg,
    hmmv2::params::PHMMParams,
    kmer::VecKmer,
    multi_dbg::MultiDbg,
    utils::{check_memory_usage, timer},
};
use std::io::Write;

/// dbgphmm experiment tools to generate synthetic genomes and reads.
///
#[derive(Parser, Debug)]
#[clap(author, about, version = env!("GIT_HASH"))]
struct Opts {
    /// k-mer size of initial DBG constructed from reads.
    /// if --datset-only, it will not affect the output.
    #[clap(short)]
    k: usize,
    /// Minimum occurrence of k-mers in reads
    /// k-mers that occurs less than min-count is removed from DBG
    #[clap(short = 'm', default_value_t = 2)]
    min_count: usize,
    /// Minimum occurrence of dead end k-mers in read
    /// All dead end k-mers whose read count is
    #[clap(short = 'M')]
    min_deadend_count: Option<usize>,
    /// Average coverage of the simulated reads.
    #[clap(short = 'C')]
    coverage: usize,
    /// Read length of the simulated reads.
    #[clap(short = 'L')]
    read_length: usize,
    /// Error rate of reads.
    /// For each bases, either mismatch/insertion/deletion will occur with a probability p.
    /// So total error rate (1 - nucleotide identity) will be 3p.
    #[clap(short = 'p')]
    p_error: f64,
    /// Random seed used to generate reads.
    #[clap(long, default_value = "0")]
    read_seed: u64,
    /// Unit length of genome
    #[clap(short = 'U')]
    unit_size: usize,
    /// How many times the unit is tandemly repeated in the genome?
    #[clap(short = 'N')]
    n_unit: usize,
    /// Length of the unique prefix/suffix of the genome.
    #[clap(short = 'E')]
    end_length: usize,
    /// inter-haplotype divergence H1
    #[clap(short = 'H')]
    hap_divergence: f64,
    /// intra-haplotype (repeat) divergence H0
    #[clap(long = "H0")]
    hap_init_divergence: f64,
    /// the number of haplotypes of the genome
    #[clap(short = 'P')]
    n_haplotypes: usize,
    /// Random seed used to generate genome (unit) sequence.
    #[clap(long, default_value = "0")]
    genome_seed: u64,
    /// Specify genome sequences as FASTA, to generate reads and DBG.
    /// If use this, the options to generate synthetic genome sequences such as -U,-N,-E,-H,--H0,-P will not affect the output.
    #[clap(long)]
    genome_fasta: Option<std::path::PathBuf>,
    /// prefix of the outputs
    #[clap(long)]
    output_prefix: std::path::PathBuf,
    /// Generate dataset (genome and reads) only, and do not construct initial DBG.
    #[clap(long)]
    dataset_only: bool,
    /// Dump raw dbg containing read count information for manual inspection.
    /// size of the dbg can be large.
    #[clap(long)]
    dump_raw_dbg: bool,
}

// -U 10000 -N 10 -E 200 -P 2 -H 0.01
// -C 20 -p 0.01 -L 20000

#[allow(unused_must_use)]
fn main() {
    let opts: Opts = Opts::parse();
    let mut log_file = std::fs::File::create(&opts.output_prefix).unwrap();

    writeln!(&mut log_file, "# started_at={}", chrono::Local::now());
    writeln!(&mut log_file, "# git_hash={}", env!("GIT_HASH"));
    writeln!(&mut log_file, "# opts={:?}", opts);
    check_memory_usage();

    let genome = if let Some(genome_fasta) = &opts.genome_fasta {
        genome::Genome::from_fasta(genome_fasta).expect("genome fasta is invalid")
    } else {
        // generate genome
        genome::tandem_repeat_polyploid_with_unique_homo_ends(
            opts.unit_size,
            opts.n_unit,
            opts.genome_seed,
            opts.hap_init_divergence,
            opts.genome_seed.wrapping_add(1),
            opts.end_length,
            opts.n_haplotypes,
            opts.hap_divergence,
            opts.genome_seed,
        )
    };
    let param = PHMMParams::uniform(opts.p_error);
    let (dataset, t) = timer(|| {
        generate_dataset(
            genome,
            opts.read_seed,
            opts.coverage,
            opts.read_length,
            ReadType::FragmentWithRevComp,
            param,
        )
    });
    writeln!(&mut log_file, "# dataset created in {}ms", t);

    // dataset.show_reads_with_genome();
    let (_, t) = timer(|| dataset.to_json_file(opts.output_prefix.with_extension("json")));
    writeln!(&mut log_file, "# dataset dumped in {}ms", t);
    dataset.to_json_file(opts.output_prefix.with_extension("json"));
    dataset.to_genome_fasta(opts.output_prefix.with_extension("genome.fa"));
    dataset.to_reads_fasta(opts.output_prefix.with_extension("reads.fa"));
    dataset.to_reads_sam(opts.output_prefix.with_extension("reads.sam"));

    if opts.dump_raw_dbg {
        let hd: HashDbg<VecKmer> = HashDbg::from_fragment_seqs(opts.k, dataset.reads());
        hd.inspect_with_genome(dataset.genome());
        hd.to_gfa_file_with_genome(
            opts.output_prefix.with_extension("rawdbg.gfa"),
            dataset.genome(),
        );
    }

    if opts.dataset_only {
        return;
    }

    let min_deadend_count = opts
        .min_deadend_count
        .unwrap_or((dataset.coverage() / 4.0) as usize);
    let (mdbg, t) = timer(|| {
        MultiDbg::create_draft_from_dataset_with(
            opts.k,
            &dataset,
            opts.min_count,
            min_deadend_count,
        )
    });
    writeln!(&mut log_file, "# draft dbg created in {}ms", t);
    mdbg.to_gfa_file(opts.output_prefix.with_extension("gfa"));
    mdbg.to_dbg_file(opts.output_prefix.with_extension("dbg"));

    let paths_true = mdbg
        .paths_from_styled_seqs(dataset.genome())
        .unwrap_or_else(|notfound| {
            writeln!(&mut log_file, "# kmer notfound {}", notfound);
            panic!("some true kmers notfound {}", notfound);
        });

    mdbg.to_paths_file(opts.output_prefix.with_extension("paths"), &paths_true);
    println!("DBG has all true kmers");

    writeln!(&mut log_file, "# finished_at={}", chrono::Local::now());
}
