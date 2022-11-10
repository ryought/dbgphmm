use clap::Clap;
use dbgphmm::e2e::{generate_dataset, ReadType};
use dbgphmm::genome;
use dbgphmm::prelude::*;

#[derive(Clap)]
struct Opts {
    #[clap(short = 'c', default_value = "20")]
    coverage: usize,
    #[clap(short = 'p', default_value = "0.003")]
    p_error: f64,
    #[clap(short, long)]
    simple: bool,
}

fn main() {
    let opts: Opts = Opts::parse();
    let (genome, genome_size) = if opts.simple {
        genome::simple(100, 5)
    } else {
        genome::tandem_repeat_haploid(20, 5, 0.01, 0, 0)
        // genome::tandem_repeat_haploid(50, 4, 0.05, 0, 0);
    };
    let coverage = opts.coverage;
    let param = PHMMParams::uniform(opts.p_error);
    let dataset = generate_dataset(
        genome,
        genome_size,
        0,
        coverage,
        2000,
        ReadType::FullLength,
        param,
    );
    let json = serde_json::to_string_pretty(&dataset).unwrap();
    println!("{}", json);
}
