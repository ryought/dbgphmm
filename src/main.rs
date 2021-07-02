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
    Sample(Sample),
    CalcProb(CalcProb),
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
}

/// Calculate probability that the model produces the reads
#[derive(Clap)]
struct CalcProb {
    /// dbg fasta file
    dbg_fa: String,
    /// reads fasta file
    reads_fa: String,
}

/// Optimize the model to fit the reads
#[derive(Clap)]
struct Optimize {
    /// dbg fasta file
    dbg_fa: String,
    /// reads fasta file
    reads_fa: String,
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
        SubCommand::Sample(t) => {
            cli::sample(t.dbg_fa, t.length, t.n_reads, k, t.seed, param);
        }
        SubCommand::CalcProb(t) => {
            cli::calc_prob(t.dbg_fa, t.reads_fa, k, param);
        }
        SubCommand::Optimize(t) => {
            cli::optimize(t.dbg_fa, t.reads_fa, k, param);
        }
        SubCommand::Sandbox(t) => {
            cli::sandbox3();
        }
    }

    // println!("emmissions: {:?}", es);
    // hmm::base::test_random();
    // hmm::testing::test_static();
}
