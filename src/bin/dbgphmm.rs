use clap::{Parser, Subcommand};
use dbgphmm::{
    common::collection::Reads,
    hashdbg::HashDbg,
    hmmv2::params::PHMMParams,
    kmer::VecKmer,
    multi_dbg::{posterior::test::test_inference_from_dbg, MultiDbg},
    prob::Prob,
};

#[derive(Parser, Debug)]
#[clap(author, about, version = env!("GIT_HASH"))]
struct Opts {
    #[clap(subcommand)]
    command: Commands,
    /// Number of threads used in parallel posterior sampling.
    /// Use all threads if not specified.
    #[clap(short = 't')]
    n_threads: Option<usize>,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Construct draft DBG from reads
    RawDbg {
        /// k of DBG
        #[clap(short = 'k')]
        k: usize,
        /// Minimum occurrence of k-mers in read
        #[clap(short = 'm', default_value_t = 2)]
        min_count: usize,
        /// Minimum occurrence of deadend k-mers in read
        #[clap(short = 'M')]
        min_deadend_count: usize,
        /// Input read FASTA filename
        read_fasta: std::path::PathBuf,
        /// Output GFA of dbg filename
        #[clap(short, long)]
        gfa_output: std::path::PathBuf,
    },
    /// Construct draft DBG from reads
    Draft {
        /// k of DBG
        #[clap(short = 'k')]
        k: usize,
        /// Minimum occurrence of k-mers in read
        #[clap(short = 'm', default_value_t = 2)]
        min_count: usize,
        /// Minimum occurrence of deadend k-mers in read
        #[clap(short = 'M')]
        min_deadend_count: usize,
        /// Expected error rate of reads.
        /// If not specified, it will use p=0.001 (0.1%) HiFi error rate.
        #[clap(short = 'p', default_value_t = 0.001)]
        p_error: f64,
        /// Expected size (total number of bases) of target genome
        #[clap(short = 'G')]
        genome_size: usize,
        /// Expected number of haplotypes in target genome if known
        #[clap(short = 'P')]
        n_haplotypes: Option<usize>,
        /// Input read FASTA filename
        read_fasta: std::path::PathBuf,
        /// Output dbg filename
        #[clap(short, long)]
        dbg_output: std::path::PathBuf,
        /// Output GFA of dbg filename
        #[clap(short, long)]
        gfa_output: Option<std::path::PathBuf>,
    },
    ///
    Infer {
        /// Filename of initial DBG
        #[clap(short, long)]
        dbg_input: std::path::PathBuf,
        /// Prefix of output files used as a working directory
        #[clap(short, long)]
        output_prefix: std::path::PathBuf,
        /// Target k of DBG
        #[clap(short = 'K')]
        k_max: usize,
        /// Expected size (total number of bases) of target genome
        #[clap(short = 'G')]
        genome_size: usize,
        /// Expected size (total number of bases) sigma (standard deviation)
        /// of target genome
        #[clap(short = 'S')]
        genome_size_sigma: usize,
        /// Input read FASTA filename
        read_fasta: std::path::PathBuf,
        /// Expected error rate of reads.
        /// If not specified, it will use p=0.001 (0.1%) HiFi error rate.
        #[clap(short = 'p', default_value_t = 0.001)]
        p_error: f64,
        /// Error rate of reads while inference
        #[clap(short = 'e', default_value_t = 0.00001)]
        p_infer: f64,
        /// If probability of copy number being zero is above p0,
        /// the k-mer will be regarded as zero-copy and purged.
        /// In default, 0.8 will be used.
        /// Larger p0, bigger the graph.
        #[clap(long, default_value_t = 0.8)]
        p0: f64,
        /// Maximum number of iteration of posterior sampling of single k
        #[clap(short = 'I', default_value = "50")]
        max_iter: usize,
        #[clap(short = 'c', default_value = "1000")]
        max_cycle_size: usize,
    },
}

fn main() {
    let opts: Opts = Opts::parse();
    println!("# started_at={}", chrono::Local::now());
    println!("# n_threads={}", rayon::current_num_threads());
    println!("# git_hash={}", env!("GIT_HASH"));
    println!("# opts={:?}", opts);

    if let Some(n_threads) = opts.n_threads {
        println!("# n_threads_specified={}", n_threads);
        rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build_global()
            .unwrap();
        println!("# n_threads={}", rayon::current_num_threads());
    }

    match &opts.command {
        Commands::RawDbg {
            k,
            min_count,
            min_deadend_count,
            read_fasta,
            gfa_output,
        } => {
            let reads = Reads::from_fasta(read_fasta).unwrap();
            let mut hd: HashDbg<VecKmer> = HashDbg::from_fragment_seqs(*k, &reads);
            hd.remove_rare_kmers(*min_count);
            hd.remove_deadends(*min_deadend_count);
            hd.to_gfa_file(gfa_output);
        }
        Commands::Draft {
            k,
            min_count,
            min_deadend_count,
            p_error,
            genome_size,
            n_haplotypes,
            read_fasta,
            dbg_output,
            gfa_output,
        } => {
            let reads = Reads::from_fasta(read_fasta).unwrap();
            println!("n_reads={}", reads.len());
            let p_error = Prob::from_prob(*p_error);
            let dbg = MultiDbg::create_draft_from_reads_v2(
                *k,
                &reads,
                p_error,
                *genome_size,
                *n_haplotypes,
                *min_count,
                *min_deadend_count,
            );

            // output
            dbg.to_dbg_file(dbg_output);
            if let Some(gfa_output) = gfa_output {
                dbg.to_gfa_file(gfa_output);
            }
        }
        Commands::Infer {
            dbg_input,
            output_prefix,
            k_max,
            genome_size,
            genome_size_sigma,
            read_fasta,
            p_error,
            p_infer,
            p0,
            max_iter,
            max_cycle_size,
        } => {
            let reads = Reads::from_fasta(read_fasta).unwrap();
            let dbg = MultiDbg::from_dbg_file(dbg_input);
            test_inference_from_dbg(
                dbg,
                &reads,
                *k_max,
                PHMMParams::uniform(*p_infer),
                PHMMParams::uniform(*p_error),
                *genome_size,
                *genome_size_sigma,
                *max_iter,
                *max_cycle_size,
                Prob::from_prob(*p0),
                output_prefix,
                None,
                None,
            );
        }
        _ => {}
    }
    println!("# finished_at={}", chrono::Local::now());
}
