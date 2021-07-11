use crate::dbg::DBG;
use crate::hmm::base::PHMM;
use crate::hmm::params::PHMMParams;
use crate::hmm::sampler::PHMMSampler;
use crate::kmer::kmer::Kmer;
use crate::prob::Prob;
use crate::*;
use clap::{AppSettings, ArgEnum, Clap};
use log::{info, warn};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::prelude::*;

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
    Compare(Compare),
    Sample(Sample),
    Forward(Forward),
    Benchmark(Benchmark),
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

/// Compare two dbg, e.g. between two refs or ref/reads
#[derive(Clap)]
struct Compare {
    /// self dbg fasta (e.g. reference fasta)
    self_dbg_fa: String,
    /// other dbg fasta (e.g. reads fasta)
    other_dbg_fa: String,
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
    /// Use rayon to parallel calculation
    #[clap(long)]
    parallel: bool,
}

#[derive(Clap)]
enum Optimizer {
    Annealer(Annealer),
    Grad(Grad),
}

#[derive(ArgEnum, Debug)]
pub enum InitStateType {
    Zero,
    True,
    ReadCount,
    Random,
}

/// Benchmark optimizers of the model to fit the reads with true dbg file
#[derive(Clap)]
struct Benchmark {
    #[clap(subcommand)]
    optimizer: Optimizer,
    /// reads fasta file
    reads_fa: String,
    /// true dbg fasta file for validation
    true_dbg_fa: String,
    /// std. dev. of genome size in prior distribution
    #[clap(short = 'V', long)]
    genome_size_std_var: u32,
    /// Only uses prior score as a state score (not forward score)
    #[clap(long)]
    prior_only: bool,
    /// initial state that optimization starts from
    #[clap(long, arg_enum, default_value = "zero")]
    init_state: InitStateType,
    /// Dump seqs as FASTA after optimizing
    #[clap(long)]
    dump_seqs: bool,
    /// Use rayon to parallel calculation
    #[clap(long)]
    parallel: bool,
    /// random seed
    #[clap(short = 's', long, default_value = "11")]
    seed: u64,
}

/// Optimize the model without answer
#[derive(Clap)]
struct Optimize {
    #[clap(subcommand)]
    optimizer: Optimizer,
    /// reads fasta file
    reads_fa: String,
    /// average of genome size in prior distribution
    #[clap(short = 'M', long)]
    genome_size_ave: u32,
    /// std. dev. of genome size in prior distribution
    #[clap(short = 'V', long)]
    genome_size_std_var: u32,
    /// Only uses prior score as a state score (not forward score)
    #[clap(long)]
    prior_only: bool,
    /// Start from true copy numbers infered from true_dbg_fa
    #[clap(long)]
    start_from_true_copy_nums: bool,
    /// Dump seqs as FASTA after optimizing
    #[clap(long)]
    dump_seqs: bool,
    /// Use rayon to parallel calculation
    #[clap(long)]
    parallel: bool,
    /// random seed
    #[clap(short = 's', long, default_value = "11")]
    seed: u64,
}

/// Optimize by simulated annealing
#[derive(Clap)]
struct Annealer {
    /// initial temperature in simulated annealing
    #[clap(short = 'T', long, default_value = "1.0")]
    init_temp: f64,
    /// cooling rate in simulated annealing
    #[clap(short = 'R', long, default_value = "0.8")]
    cooling_rate: f64,
    /// cooling rate in simulated annealing
    #[clap(short = 'I', long, default_value = "100")]
    n_iteration: u64,
}

/// Optimize by gradient
#[derive(Clap)]
struct Grad {
    /// max iteration number
    #[clap(short = 'I', long, default_value = "10")]
    max_iteration: u64,
}

/// Sandbox for debugging
#[derive(Clap)]
struct Sandbox {}

pub fn main() {
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
        SubCommand::Generate(o) => {
            cli::generate(o);
        }
        SubCommand::Stat(o) => {
            cli::stat(o, k);
        }
        SubCommand::Compare(o) => {
            cli::compare(o, k);
        }
        SubCommand::Sample(o) => {
            cli::sample(o, k, param);
        }
        SubCommand::Forward(o) => {
            cli::forward(o, k, param);
        }
        SubCommand::Benchmark(o) => {
            cli::benchmark(o, k, param);
        }
        SubCommand::Optimize(o) => {
            cli::optimize(o, k, param);
        }
        SubCommand::Sandbox(o) => {
            cli::sandbox(o);
        }
    }
}

fn generate(opts: Generate) {
    let v = random_seq::generate(opts.length, opts.seed);
    println!(">randseq");
    println!("{}", std::str::from_utf8(&v).unwrap());
}

fn stat(opts: Stat, k: usize) {
    let seqs = io::fasta::parse_seqs(&opts.dbg_fa);
    let d = dbg::DbgHash::from_seqs(&seqs, k);
    info!("{:?}", d.as_degree_stats());
    let cdbg = compressed_dbg::CompressedDBG::from(&d, k);
    info!("{:?}", cdbg.as_dbg_stats());
    info!("{:?}", cdbg.as_degree_stats());
    let copy_nums_true: Vec<u32> = cdbg
        .iter_nodes()
        .map(|v| {
            let kmer = cdbg.kmer(&v);
            d.find(kmer)
        })
        .collect();

    info!("{:?}", cdbg.as_copy_num_stats(&copy_nums_true));
    // let cn_zero: Vec<u32> = vec![0; cdbg.n_kmers()];
    // info!("{:?}", cdbg.as_copy_num_stats(&cn_zero));
    // info!("{:?}", copy_nums_true);
    info!("cycle {}", cdbg.as_cycle_histogram());
    info!("{:?}", cdbg.as_cycle_summary_stats());

    for i in 0..cdbg.n_cycles() {
        info!("{:?}", cdbg.as_cycle_stats(i));
    }

    let all_stats = cdbg.as_all_stats(&copy_nums_true);
    let json = serde_json::to_string_pretty(&all_stats).unwrap();
    // info!("{:#?}", cdbg.as_all_stats(&copy_nums_true));
    println!("{}", json);
}

fn compare(opts: Compare, k: usize) {
    let self_seqs = io::fasta::parse_seqs(&opts.self_dbg_fa);
    let self_dbg = dbg::DbgHash::from_seqs(&self_seqs, k);

    let other_seqs = io::fasta::parse_seqs(&opts.other_dbg_fa);
    let other_dbg = dbg::DbgHash::from_seqs(&other_seqs, k);

    let result = self_dbg.compare_dbg(&other_dbg);
    let json = serde_json::to_string_pretty(&result).unwrap();
    println!("{}", json);
}

fn sample(opts: Sample, k: usize, param: PHMMParams) {
    let seqs = io::fasta::parse_seqs(&opts.dbg_fa);
    let (cdbg, copy_nums) = compressed_dbg::CompressedDBG::from_seqs(&seqs, k);
    let phmm = hmm::cdbg::CDbgPHMM::new(&cdbg, copy_nums);
    info!("{:?}", cdbg.as_degree_stats());
    let from = if opts.start_from_head {
        let head = cdbg
            .heads()
            .first()
            .unwrap_or_else(|| panic!("Cannot find head node"))
            .clone();
        Some(head)
    } else {
        None
    };
    info!("from={:?}", from);

    // println!("{}", d.as_dot());
    // println!("{}", d.dbg.as_dot());
    // let mut rng = Xoshiro256PlusPlus::seed_from_u64(seed);
    let mut infos: Vec<hmm::sampler::SampleInfo> = Vec::new();
    for i in 0..opts.n_reads {
        // let seed_for_a_read: u64 = rng.gen();
        let seed_for_a_read = opts.seed + i as u64;
        let (seq, info) = phmm.sample(&param, opts.length, seed_for_a_read, from);
        // output fasta
        let id = format!("{},{}", opts.length, seed_for_a_read);
        io::fasta::dump_seq(&id, &seq, Some(&info.to_string()));
        // store info in vec
        infos.push(info);
    }
    info!("{:?}", hmm::sampler::sum_sample_infos(&infos));
}

fn forward(opts: Forward, k: usize, param: PHMMParams) {
    let seqs = io::fasta::parse_seqs(&opts.dbg_fa);
    let (cdbg, copy_nums) = compressed_dbg::CompressedDBG::from_seqs(&seqs, k);
    let phmm = hmm::cdbg::CDbgPHMM::new(&cdbg, copy_nums);

    let reads = io::fasta::parse_seqs(&opts.reads_fa);

    if opts.parallel {
        let p_total: prob::Prob = reads
            .par_iter()
            .map(|read| phmm.forward_prob(&param, read))
            .product();
        println!("#total\t{}", p_total.to_log_value());
    } else {
        let mut ps: Vec<prob::Prob> = Vec::new();
        for (i, read) in reads.iter().enumerate() {
            let p = phmm.forward_prob(&param, read);
            println!("{}\t{}", i, p.to_log_value());
            ps.push(p);
            // let p = phmm.backward_prob(&param, read);
            // println!("backward prob : {}", p);
        }
        let p_total: prob::Prob = ps.iter().product();
        println!("#total\t{}", p_total.to_log_value());
    }
}

/// Experiments of optimizer
/// 1. construct dbg from reads
/// 2. determine (true) copy_nums from fa
/// 3. optimize
fn benchmark(opts: Benchmark, k: usize, param: PHMMParams) {
    let reads = io::fasta::parse_seqs(&opts.reads_fa);
    let (cdbg, copy_nums_read) = compressed_dbg::CompressedDBG::from_seqs(&reads, k);

    let seqs = io::fasta::parse_seqs(&opts.true_dbg_fa);
    let copy_nums_true = cdbg
        .true_copy_nums_from_seqs(&seqs, k)
        .unwrap_or_else(|| panic!("True copy_nums is not in read cdbg"));
    let true_size = cdbg.total_emitable_copy_num(&copy_nums_true);
    info!("true_size={}", true_size);

    let true_state = optimizer::cdbg::CDbgState::new(
        &cdbg,
        copy_nums_true.clone(),
        cdbg.cycle_vec_from_copy_nums(&copy_nums_true),
        true_size,
        opts.genome_size_std_var,
        if opts.prior_only { None } else { Some(&reads) },
        param.clone(),
        opts.parallel,
    );

    match opts.optimizer {
        Optimizer::Annealer(opts_annealer) => {
            let mut rng = Xoshiro256PlusPlus::seed_from_u64(opts.seed);
            let a = optimizer::annealer::Annealer::new(
                opts_annealer.init_temp,
                opts_annealer.cooling_rate,
            );

            match opts.init_state {
                InitStateType::Zero => {
                    // test run from true
                    a.run_with_log(&mut rng, true_state, 1);
                    // real run from zero
                    let init_state = optimizer::cdbg::CDbgState::init(
                        &cdbg,
                        true_size,
                        opts.genome_size_std_var,
                        if opts.prior_only { None } else { Some(&reads) },
                        param.clone(),
                        opts.parallel,
                    );
                    let history = a.run_with_log(&mut rng, init_state, opts_annealer.n_iteration);
                    // println!("{:?}", history.last().unwrap().copy_nums);
                    // println!("{:?}", copy_nums_true);
                    let copy_nums_final = &history.last().unwrap().copy_nums;
                    if opts.dump_seqs {
                        for (i, seq) in cdbg.to_seqs(copy_nums_final).iter().enumerate() {
                            let id = format!("{}", i);
                            io::fasta::dump_seq(&id, &seq, None);
                        }
                    }
                }
                InitStateType::True => {
                    // real run from true
                    let history = a.run_with_log(&mut rng, true_state, opts_annealer.n_iteration);
                    let copy_nums_final = &history.last().unwrap().copy_nums;
                    if opts.dump_seqs {
                        for (i, seq) in cdbg.to_seqs(copy_nums_final).iter().enumerate() {
                            let id = format!("{}", i);
                            io::fasta::dump_seq(&id, &seq, None);
                        }
                    }
                }
                _ => panic!("not implemented"),
            }

            if opts.dump_seqs {
                for (i, seq) in cdbg.to_seqs(&copy_nums_true).iter().enumerate() {
                    let id = format!("t{}", i);
                    io::fasta::dump_seq(&id, &seq, None);
                }
            }
        }
        Optimizer::Grad(opts_grad) => {
            let g = optimizer::grad::GradientDescent::new(opts_grad.max_iteration, true);
            match opts.init_state {
                InitStateType::True => {
                    let history = g.run(true_state);
                }
                InitStateType::Zero => {
                    let init_state = optimizer::cdbg::CDbgState::init(
                        &cdbg,
                        true_size,
                        opts.genome_size_std_var,
                        if opts.prior_only { None } else { Some(&reads) },
                        param.clone(),
                        opts.parallel,
                    );
                    let history = g.run(init_state);
                }
                _ => panic!("not implemented"),
            }
        }
    }
}

fn optimize(opts: Optimize, k: usize, param: PHMMParams) {
    println!("not implemented!");
}

fn sandbox(opts: Sandbox) {
    let x: Vec<(i32, f32)> = vec![(1, 0.7), (2, 0.1), (3, 0.9), (4, f32::NAN)];
    println!("{:?}", x);
    let max = x
        .into_iter()
        .filter(|(_, x)| !x.is_nan())
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
        .unwrap();
    println!("{:?}", max);

    // stats::test();
    /*
    let mut reads: Vec<Vec<u8>> = Vec::new();
    // reads.push(b"ATCGATTCGATCGATTCGATAGATCG".to_vec());
    reads.push(b"AGGCTAGTAAAAAAAAAAAAAATCGATCTTTCGATCG".to_vec());
    reads.push(b"GGATAGTTCGATCTG".to_vec());
    reads.push(b"GGCTAGTTCGATCGG".to_vec());
    let (cdbg, copy_nums) = compressed_dbg::CompressedDBG::from_seqs(&reads, 8);
    // println!("{}", cdbg.as_dot_with_copy_nums(&copy_nums));
    for (i, seq) in cdbg.to_seqs(&copy_nums).iter().enumerate() {
        let id = format!("{}", i);
        io::fasta::dump_seq(&id, &seq, None);
    }
    */
}
