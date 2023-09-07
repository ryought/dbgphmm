use clap::{Parser, Subcommand};
use dbgphmm::{
    common::collection::{Reads, Seq},
    dbg::draft::EndNodeInference,
    distribution::kmer_coverage,
    hashdbg::HashDbg,
    kmer::VecKmer,
    multi_dbg::MultiDbg,
    prob::Prob,
};

#[derive(Parser, Debug)]
#[clap(author, about, version = env!("GIT_HASH"))]
struct Opts {
    #[clap(subcommand)]
    command: Commands,
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
        /// Expected number of haplotypes in target genome if known
        #[clap(short = 'P')]
        n_haplotypes: Option<usize>,
        /// Input read FASTA filename
        read_fasta: std::path::PathBuf,
    },
}

fn main() {
    let opts: Opts = Opts::parse();
    println!("# started_at={}", chrono::Local::now());
    println!("# n_threads={}", rayon::current_num_threads());
    println!("# git_hash={}", env!("GIT_HASH"));
    println!("# opts={:?}", opts);
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
        _ => {}
    }
    println!("# finished_at={}", chrono::Local::now());
}
