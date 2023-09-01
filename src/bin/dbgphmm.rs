use clap::{Parser, Subcommand};

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
        dbg: std::path::PathBuf,
    },
    ///
    Infer {
        /// Filename of initial DBG
        #[clap(short, long)]
        dbg: std::path::PathBuf,
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
            dbg,
        } => {
            println!("draft..");
        }
        _ => {}
    }
    println!("# finished_at={}", chrono::Local::now());
}
