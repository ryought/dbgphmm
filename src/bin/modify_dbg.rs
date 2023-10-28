use clap::Parser;
use dbgphmm::multi_dbg::MultiDbg;

#[derive(Parser, Debug)]
#[clap(author, about, version = env!("GIT_HASH"))]
struct Opts {
    #[clap(long)]
    dbg_in: std::path::PathBuf,
    #[clap(long)]
    inspect: std::path::PathBuf,
    #[clap(long)]
    sample_id: usize,
    #[clap(long)]
    dbg_out: std::path::PathBuf,
}

fn main() {
    let opts: Opts = Opts::parse();
    println!("# started_at={}", chrono::Local::now());
    println!("# n_threads={}", rayon::current_num_threads());
    println!("# git_hash={}", env!("GIT_HASH"));
    println!("# opts={:?}", opts);
    println!("loading dbg");
    let mut dbg = MultiDbg::from_dbg_file(&opts.dbg_in);
    println!("loading inspect");
    let posterior = dbg.from_inspect_file(&opts.inspect).unwrap();
    // for sample in posterior.samples {
    //     println!("{:?}", sample);
    // }
    dbg.set_copy_nums(&posterior.samples[opts.sample_id].copy_nums);
    dbg.to_dbg_file(&opts.dbg_out).unwrap();
}
