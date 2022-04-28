//!
//! mock genome generation functions
//!
//! * `simple`
//! * `simple_diploid`
//! * `tandem_repeat`
//! * `tandem_repeat_diploid`
//!
use crate::common::{sequence_to_string, Genome, Reads, Seq, Sequence};
use crate::random_seq::{generate, random_mutation, tandem_repeat, MutationProfile};

/// tandem repeat haploid genome
///
pub fn generate_tandem_repeat_haploid(
    unit_size: usize,
    n_unit: usize,
    divergence_init: f64,
    unit_seed: u64,
    hap_seed: u64,
) -> (Genome, usize) {
    let genome_size = unit_size * n_unit;
    let unit = generate(unit_size, unit_seed);
    let tandem_repeat = tandem_repeat(&unit, n_unit);
    let (hap_a, _) = random_mutation(
        &tandem_repeat,
        MutationProfile::uniform(divergence_init),
        hap_seed,
    );
    println!("{}", tandem_repeat.to_str());
    println!("{}", hap_a.to_str());
    (vec![hap_a], genome_size)
}

// fn generate_diploid_from_haploid()

pub fn generate_tandem_repeat_diploid(
    unit_size: usize,
    n_unit: usize,
    divergence_init: f64,
    unit_seed: u64,
    hap_seed: u64,
    divergence_between_haplotypes: f64,
    div_seed: u64,
) -> (Genome, usize) {
    let (mut hap, hap_genome_size) =
        generate_tandem_repeat_haploid(unit_size, n_unit, divergence_init, unit_seed, hap_seed);
    let hap_a = hap.remove(0);
    let (hap_b, _) = random_mutation(
        &hap_a,
        MutationProfile::uniform(divergence_between_haplotypes),
        div_seed,
    );
    let genome_size = hap_genome_size * 2;
    println!("{}", hap_a.to_str());
    println!("{}", hap_b.to_str());
    (vec![hap_a, hap_b], genome_size)
}
