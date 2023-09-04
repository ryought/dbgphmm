use clap::{Parser, Subcommand};
use dbgphmm::{
    common::collection::{Reads, Seq},
    dbg::draft::EndNodeInference,
    hashdbg::HashDbg,
    kmer::VecKmer,
    multi_dbg::MultiDbg,
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
    Draft {
        /// k of DBG
        #[clap(short = 'k')]
        k: usize,
        /// Minimum occurrence of k-mers in read
        #[clap(short = 'm', default_value_t = 2)]
        min_count: usize,
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
        Commands::Draft {
            k,
            min_count,
            p_error,
            genome_size,
            n_haplotypes,
            read_fasta,
            dbg_output,
            gfa_output,
        } => {
            println!("draft..");
            let reads = Reads::from_fasta(read_fasta).unwrap();
            println!("n_reads={}", reads.len());
            let mut hd: HashDbg<VecKmer> = HashDbg::from_fragment_seqs(*k, &reads);
            hd.remove_rare_kmers(*min_count);
            if let Some(gfa_output) = gfa_output {
                hd.to_gfa_file(gfa_output);
            }

            /*
            let dbg = MultiDbg::create_draft_from_reads(
                *k,
                &reads,
                reads.coverage(*genome_size),
                reads.average_length(),
                *p_error,
                &EndNodeInference::Auto,
            );
            dbg.to_dbg_file(dbg_output);
            if let Some(gfa_output) = gfa_output {
                dbg.to_gfa_file(gfa_output);
            }
            */
        }
        _ => {}
    }
    println!("# finished_at={}", chrono::Local::now());
}
