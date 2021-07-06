use clap::{AppSettings, Clap};
use dbgphmm::hmm::base::PHMM;
use dbgphmm::prob::Prob;
use dbgphmm::*;
use log::{info, warn};
use std::io::prelude::*;

/// de bruijn graph + profile HMM optimization package
#[derive(Clap)]
#[clap(version = "0.1", author = "ryought <ryonakabayashi@gmail.com>")]
#[clap(setting = AppSettings::ColoredHelp)]
struct Opts {
    #[clap(subcommand)]
    subcmd: SubCommand,
    /// kmer size
    #[clap(short = 'k', long, default_value = "8")]
    kmer_size: usize,
    /// mismatch probability
    #[clap(long, default_value = "0.01")]
    p_mismatch: f64,
    /// gap open probability
    #[clap(long, default_value = "0.01")]
    p_gap_open: f64,
    /// gap extension probability
    #[clap(long, default_value = "0.01")]
    p_gap_ext: f64,
    /// end probability
    #[clap(long, default_value = "0.00001")]
    p_end: f64,
    /// number of consecutive dels to consider
    #[clap(long, default_value = "3")]
    n_max_gaps: u32,
    /// Print debug info
    #[clap(short, long)]
    debug: bool,
}

#[derive(Clap)]
enum SubCommand {
    Generate(Generate),
    Stat(Stat),
    ReadStat(ReadStat),
    Sample(Sample),
    Forward(Forward),
    Optimize(Optimize),
    Sandbox(Sandbox),
}

/// Generate a random sequence
#[derive(Clap)]
struct Generate {
    /// sequence length to generate
    #[clap(short = 'l', long)]
    length: usize,
    /// random seed
    #[clap(short = 's', long, default_value = "0")]
    seed: u64,
}

/// Show de bruijn graph statistics
#[derive(Clap)]
struct Stat {
    /// dbg fasta file
    dbg_fa: String,
}

/// Statistics of reads
#[derive(Clap)]
struct ReadStat {
    /// dbg fasta file
    dbg_fa: String,
    /// read fasta file
    reads_fa: String,
}

/// Sample reads from the model
#[derive(Clap)]
struct Sample {
    /// dbg fasta file
    dbg_fa: String,
    /// read length to generate
    #[clap(short = 'l', long)]
    length: u32,
    /// number of reads
    #[clap(short = 'n', long)]
    n_reads: u32,
    /// random seed
    #[clap(short = 's', long, default_value = "0")]
    seed: u64,
    /// Require all reads to start from the head node
    #[clap(long)]
    start_from_head: bool,
}

/// Calculate probability that the model produces the reads
#[derive(Clap)]
struct Forward {
    /// dbg fasta file
    dbg_fa: String,
    /// reads fasta file
    reads_fa: String,
}

/// Optimize the model to fit the reads
#[derive(Clap)]
struct Optimize {
    /// reads fasta file
    reads_fa: String,
    /// true dbg fasta file for validation
    #[clap(long)]
    true_dbg_fa: Option<String>,
    /// initial temperature in simulated annealing
    #[clap(short = 'T', long, default_value = "1.0")]
    init_temp: f64,
    /// cooling rate in simulated annealing
    #[clap(short = 'R', long, default_value = "0.8")]
    cooling_rate: f64,
    /// cooling rate in simulated annealing
    #[clap(short = 'I', long, default_value = "100")]
    n_iteration: u64,
    /// average of genome size in prior distribution
    #[clap(short = 'M', long)]
    genome_size_ave: u32,
    /// std. dev. of genome size in prior distribution
    #[clap(short = 'V', long)]
    genome_size_std_var: u32,
    /// Only uses prior score as a state score (not forward score)
    #[clap(long)]
    prior_only: bool,
}

/// Sandbox for debugging
#[derive(Clap)]
struct Sandbox {}

fn main() {
    // enable logger
    env_logger::init();

    // parse options
    let opts: Opts = Opts::parse();
    let param = hmm::params::PHMMParams::new(
        prob::Prob::from_prob(opts.p_mismatch),
        prob::Prob::from_prob(opts.p_gap_open),
        prob::Prob::from_prob(opts.p_gap_ext),
        prob::Prob::from_prob(opts.p_end),
        opts.n_max_gaps,
    );
    let k = opts.kmer_size;

    match opts.subcmd {
        SubCommand::Generate(t) => {
            cli::generate(t.length, t.seed);
        }
        SubCommand::Stat(t) => {
            cli::stat(t.dbg_fa, k);
        }
        SubCommand::ReadStat(t) => {
            cli::readstat(t.dbg_fa, t.reads_fa, k);
        }
        SubCommand::Sample(t) => {
            cli::sample(
                t.dbg_fa,
                t.length,
                t.n_reads,
                k,
                t.seed,
                param,
                t.start_from_head,
            );
        }
        SubCommand::Forward(t) => {
            cli::forward(t.dbg_fa, t.reads_fa, k, param);
        }
        SubCommand::Optimize(t) => match t.true_dbg_fa {
            Some(dbg_fa) => cli::optimize_with_answer(
                dbg_fa,
                t.reads_fa,
                k,
                param,
                t.init_temp,
                t.cooling_rate,
                t.n_iteration,
                t.genome_size_ave,
                t.genome_size_std_var,
                t.prior_only,
            ),
            None => cli::optimize(t.reads_fa, k, param),
        },
        SubCommand::Sandbox(t) => {
            cli::sandbox();
        }
    }

    // println!("emmissions: {:?}", es);
    // hmm::base::test_random();
    // hmm::testing::test_static();
}
