use bio::alignment::distance::levenshtein;
use bio::alignment::pairwise::*;
use clap::Parser;
use dbgphmm::{
    e2e::Dataset,
    multi_dbg::{MultiDbg, NeighborConfig},
};

#[derive(Parser, Debug)]
#[clap(author, about, version = env!("GIT_HASH"))]
struct Opts {
    #[clap(long)]
    dbg: std::path::PathBuf,
    #[clap(long)]
    dataset_json: std::path::PathBuf,
    #[clap(long)]
    show_alignment: bool,
}

fn main() {
    let opts: Opts = Opts::parse();
    println!("# started_at={}", chrono::Local::now());
    println!("# git_hash={}", env!("GIT_HASH"));
    println!("# opts={:?}", opts);

    let dataset = Dataset::from_json_file(opts.dataset_json);
    let dbg = MultiDbg::from_dbg_file(&opts.dbg);
    dbg.to_fasta(opts.dbg.with_extension("fa"));

    let s = dbg.to_styled_seqs();

    for (i, s) in dbg.to_styled_seqs().into_iter().enumerate() {
        for (j, g) in dataset.genome().into_iter().enumerate() {
            let x = s.seq();
            let y = g.seq();

            let mut aligner =
                Aligner::with_capacity(
                    x.len(),
                    y.len(),
                    -5,
                    -1,
                    |a: u8, b: u8| if a == b { 1 } else { -1 },
                );
            let alignment = aligner.global(x, y);

            println!("p[{}] g[{}] {}", i, j, levenshtein(x, y));

            if opts.show_alignment {
                println!("{}", alignment.pretty(x, y, 1000));
            }
        }
    }

    // TODO
    // solve bipartite graph matching and determine the best correspondence between haplotype and
    // path
}
