//!
//! mock genome generation functions
//!
//! * `simple`
//! * `simple_diploid`
//! * `tandem_repeat`
//! * `tandem_repeat_diploid`
//!
use crate::common::{sequence_to_string, Genome, Reads, Seq, Sequence, StyledSequence};
use crate::random_seq::{
    generate, join, random_mutation, random_mutation_with_rng, tandem_repeat, MutationProfile,
};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;

///
/// simple random haploid genome
///
/// ## Parameters
///
/// * genome_size: size of the genome
/// * seed: seed of the sequence generation
///
pub fn simple(genome_size: usize, seed: u64) -> (Genome, usize) {
    let genome = vec![StyledSequence::linear(generate(genome_size, seed))];
    (genome, genome_size)
}

///
/// simple random haploid circular genome
///
pub fn simple_circular(genome_size: usize, seed: u64) -> (Genome, usize) {
    let genome = vec![StyledSequence::circular(generate(genome_size, seed))];
    (genome, genome_size)
}

///
/// simple diploid genome, with two SNVs manually added.
///
/// will be deprecated
///
pub fn simple_diploid() -> (Genome, usize) {
    let genome_size_hap = 100;
    let haplotype1 = generate(100, 0);
    let mut haplotype2 = haplotype1.clone();
    haplotype2[30] = b'C';
    haplotype2[80] = b'T';
    let genome = vec![
        StyledSequence::linear(haplotype1),
        StyledSequence::linear(haplotype2),
    ];
    let genome_size = genome_size_hap * 2;
    (genome, genome_size)
}

///
/// simple diploid genome, with random SNVs
///
/// ## Parameters
///
/// * hap_size: size of the haplotype genome
/// * hap_seed: seed of haplotype sequence generation
/// * div_rate: rate of random mutation per base
/// * div_seed: seed of random mutation
///
pub fn diploid(hap_size: usize, hap_seed: u64, div_rate: f64, div_seed: u64) -> (Genome, usize) {
    let hap_a = generate(hap_size, hap_seed);
    let (hap_b, ops) = random_mutation(&hap_a, MutationProfile::uniform(div_rate), div_seed);
    // println!("ops={:?}", ops);
    let genome_size = hap_size * 2;
    (
        vec![StyledSequence::linear(hap_a), StyledSequence::linear(hap_b)],
        genome_size,
    )
}

/// tandem repeat haploid genome
///
/// -->-->-->-->-->-->-->
///  ^
///  |
/// -->-->-->-->-->-->-->
///  ^
///  |
/// -->
///
pub fn tandem_repeat_haploid(
    unit_size: usize,
    n_unit: usize,
    divergence_init: f64,
    unit_seed: u64,
    hap_seed: u64,
) -> (Genome, usize) {
    tandem_repeat_haploid_with_unique_ends(
        unit_size,
        n_unit,
        divergence_init,
        unit_seed,
        hap_seed,
        0,
    )
}

pub fn tandem_repeat_haploid_with_unique_ends(
    unit_size: usize,
    n_unit: usize,
    divergence_init: f64,
    unit_seed: u64,
    hap_seed: u64,
    end_length: usize,
) -> (Genome, usize) {
    let genome_size = unit_size * n_unit + end_length * 2;
    let unit = generate(unit_size, unit_seed);
    let tandem_repeat = tandem_repeat(&unit, n_unit);
    let (hap_a, _) = random_mutation(
        &tandem_repeat,
        MutationProfile::uniform(divergence_init),
        hap_seed,
    );
    let prefix = generate(end_length, unit_seed.wrapping_add(1));
    let suffix = generate(end_length, unit_seed.wrapping_sub(1));
    let hap_a = join(prefix, join(hap_a, suffix));
    (vec![StyledSequence::linear(hap_a)], genome_size)
}

pub fn tandem_repeat_diploid(
    unit_size: usize,
    n_unit: usize,
    divergence_init: f64,
    unit_seed: u64,
    hap_seed: u64,
    divergence_between_haplotypes: f64,
    div_seed: u64,
) -> (Genome, usize) {
    let (mut hap, hap_genome_size) =
        tandem_repeat_haploid(unit_size, n_unit, divergence_init, unit_seed, hap_seed);
    let hap_a = hap.remove(0);
    let (hap_b_seq, ops) = random_mutation(
        &hap_a.seq(),
        MutationProfile::uniform(divergence_between_haplotypes),
        div_seed,
    );
    let hap_b = StyledSequence::linear(hap_b_seq);
    let genome_size = hap_genome_size * 2;
    (vec![hap_a, hap_b], genome_size)
}

pub fn tandem_repeat_polyploid_with_unique_ends(
    unit_size: usize,
    n_unit: usize,
    divergence_init: f64,
    unit_seed: u64,
    hap_seed: u64,
    end_length: usize,
    n_haplotypes: usize,
    divergence_between_haplotypes: f64,
    div_seed: u64,
) -> (Genome, usize) {
    let (mut hap, hap_genome_size) = tandem_repeat_haploid_with_unique_ends(
        unit_size,
        n_unit,
        divergence_init,
        unit_seed,
        hap_seed,
        end_length,
    );
    let mut genome = Vec::new();
    let mut genome_size = 0;
    let hap_a = hap.remove(0);
    genome_size += hap_a.len();
    // assert_eq!(hap_a.len(), hap_genome_size);
    genome.push(hap_a.clone());

    let mut rng = Xoshiro256PlusPlus::seed_from_u64(div_seed);
    for _ in 1..n_haplotypes {
        let (hap_b_seq, ops) = random_mutation_with_rng(
            &hap_a.seq(),
            MutationProfile::uniform(divergence_between_haplotypes),
            &mut rng,
        );
        let hap_b = StyledSequence::linear(hap_b_seq);
        genome_size += hap_b.len();
        genome.push(hap_b);
    }
    (genome, genome_size)
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

    fn show_genome(genome: &Genome, genome_size: usize) {
        println!("size: {}", genome_size);
        for (i, hap) in genome.iter().enumerate() {
            println!("hap{}: {}", i, hap.to_str());
        }
    }

    #[test]
    fn e2e_genome_generation() {
        let (g, gs) = simple(100, 0);
        show_genome(&g, gs);
        assert_eq!(gs, 100);
        assert_eq!(g.len(), 1);
        assert_eq!(g, vec![StyledSequence::linear(b"CCAATTCACAAAAACCACACCTTGGCCAAGGTATCGTATCTTGTTGTTGTATGTGAAAGGGGCCCTAAGATCTGTAGCCACCATGGCTAGGGTCAAATCT".to_vec())]);

        let (g, gs) = simple(50, 5);
        show_genome(&g, gs);
        assert_eq!(gs, 50);
        assert_eq!(g.len(), 1);
        assert_eq!(
            g,
            vec![StyledSequence::linear(
                b"CGAAGATGAGAAACCCGGGAGTCGATATATTCAAACAAACGGGCGCTCCT".to_vec()
            )]
        );

        let (g, gs) = simple_diploid();
        show_genome(&g, gs);
        assert_eq!(gs, 200);
        assert_eq!(g.len(), 2);
        assert_eq!(
            g,
            vec![
                StyledSequence::linear(b"CCAATTCACAAAAACCACACCTTGGCCAAGGTATCGTATCTTGTTGTTGTATGTGAAAGGGGCCCTAAGATCTGTAGCCACCATGGCTAGGGTCAAATCT".to_vec()),
                StyledSequence::linear(b"CCAATTCACAAAAACCACACCTTGGCCAAGCTATCGTATCTTGTTGTTGTATGTGAAAGGGGCCCTAAGATCTGTAGCCATCATGGCTAGGGTCAAATCT".to_vec()),
            ]
        );

        let (g, gs) = diploid(50, 5, 0.1, 0);
        show_genome(&g, gs);
        assert_eq!(gs, 100);
        assert_eq!(g.len(), 2);
        assert_eq!(
            g,
            vec![
                StyledSequence::linear(
                    b"CGAAGATGAGAAACCCGGGAGTCGATATATTCAAACAAACGGGCGCTCCT".to_vec()
                ),
                StyledSequence::linear(
                    b"CGAAGAATGAGAAACACCGGAGTCGTATATTCCAAACAAACGGGCGCTCCT".to_vec()
                )
            ]
        );

        let (g, gs) = tandem_repeat_haploid(10, 5, 0.0, 0, 0);
        show_genome(&g, gs);
        assert_eq!(gs, 50);
        assert_eq!(g.len(), 1);
        assert_eq!(
            g,
            vec![StyledSequence::linear(
                b"CCAATTCACACCAATTCACACCAATTCACACCAATTCACACCAATTCACA".to_vec()
            ),]
        );

        let (g, gs) = tandem_repeat_haploid(10, 5, 0.1, 55, 33);
        show_genome(&g, gs);
        assert_eq!(gs, 50);
        assert_eq!(g.len(), 1);
        assert_eq!(
            g,
            vec![StyledSequence::linear(
                b"GTAAATGCGGGTAAATTGCGGGTAAATGCGGCGTAAATGCGGGGAAATCGGG".to_vec()
            )]
        );

        let (g, gs) = tandem_repeat_diploid(10, 5, 0.1, 55, 33, 0.05, 22);
        show_genome(&g, gs);
        assert_eq!(gs, 100);
        assert_eq!(g.len(), 2);
        assert_eq!(
            g,
            vec![
                StyledSequence::linear(
                    b"GTAAATGCGGGTAAATTGCGGGTAAATGCGGCGTAAATGCGGGGAAATCGGG".to_vec()
                ),
                StyledSequence::linear(
                    b"GTAAATGCGGGTAACTCTGCGGGTAAATGCGGCGTAAATGCGGGGACAATCGGG".to_vec()
                ),
            ]
        );

        let (g, gs) = tandem_repeat_haploid_with_unique_ends(10, 5, 0.0, 0, 0, 10);
        show_genome(&g, gs);
        assert_eq!(gs, 70);
        assert_eq!(g.len(), 1);
        assert_eq!(
            g,
            vec![StyledSequence::linear(
                b"TAGGACAAGCCCAATTCACACCAATTCACACCAATTCACACCAATTCACACCAATTCACACCTCACCTCA".to_vec()
            ),]
        );
    }
    #[test]
    fn genome_tandem_repeat_polyploid() {
        // ploidy=1
        let (g, gs) = tandem_repeat_polyploid_with_unique_ends(10, 5, 0.0, 0, 0, 0, 1, 0.01, 0);
        show_genome(&g, gs);
        println!("{}", gs);
        assert_eq!(
            g,
            vec![StyledSequence::linear(
                b"CCAATTCACACCAATTCACACCAATTCACACCAATTCACACCAATTCACA".to_vec()
            ),]
        );

        // ploidy=4
        let (g, gs) = tandem_repeat_polyploid_with_unique_ends(10, 5, 0.0, 0, 0, 0, 4, 0.01, 0);
        show_genome(&g, gs);
        println!("{}", gs);
        assert_eq!(
            g,
            vec![
                StyledSequence::linear(
                    b"CCAATTCACACCAATTCACACCAATTCACACCAATTCACACCAATTCACA".to_vec()
                ),
                StyledSequence::linear(
                    b"CCAATTCACACCAATTGACACCAATTCACACCAATTCACACCAATTCACA".to_vec()
                ),
                StyledSequence::linear(
                    b"CCAATTCACACCAATTCACACCAATCACACCAATTCACACCAATTCACA".to_vec()
                ),
                StyledSequence::linear(
                    b"CCAATTCACACCAAATCACACCAATTCACACCAATTCACACCAATTCACA".to_vec()
                ),
            ]
        );
    }
}
