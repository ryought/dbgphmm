use clap::Clap;
use dbgphmm::e2e::{generate_dataset, ReadType};
use dbgphmm::genome;
use dbgphmm::prelude::*;

#[derive(Clap)]
struct Opts {
    #[clap(short = 'c', default_value = "30")]
    coverage: usize,
    #[clap(short = 'p', default_value = "0.001")]
    p_error: f64,
    #[clap(short, long)]
    simple: bool,
}

fn main() {
    let opts: Opts = Opts::parse();
    let (genome, genome_size) = if opts.simple {
        // 100kb unique genome
        genome::simple(100_000, 5)
    } else {
        // 100kb genome (10k unit x 10)
        genome::tandem_repeat_haploid(10_000, 10, 0.01, 0, 0)
        // genome::tandem_repeat_haploid(50, 4, 0.05, 0, 0);
    };
    let coverage = opts.coverage;
    let param = PHMMParams::uniform(opts.p_error);
    let dataset = generate_dataset(
        genome,
        genome_size,
        0,
        coverage,
        genome_size * 2,
        ReadType::FullLength,
        param,
    );
    let json = serde_json::to_string_pretty(&dataset).unwrap();
    println!("{}", json);
}
