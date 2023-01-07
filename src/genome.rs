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

pub fn tandem_repeat_diploid_example_mut() -> (Genome, usize) {
    let (mut genome, genome_size) =
        tandem_repeat_polyploid_with_unique_ends(50, 20, 0.00, 0, 0, 50, 2, 0.00, 0);
    genome[1].seq[75] = b'A';
    (genome, genome_size)
}

pub fn tandem_repeat_diploid_example_ins() -> (Genome, usize) {
    let (mut genome, genome_size) =
        tandem_repeat_polyploid_with_unique_ends(50, 20, 0.00, 0, 0, 50, 2, 0.00, 0);
    genome[1].seq.insert(75, b'C');
    (genome, genome_size)
}

pub fn tandem_repeat_500bp() -> (Genome, usize) {
    let seed = 1;
    tandem_repeat_polyploid_with_unique_ends(10, 50, 0.0, seed, seed, 50, 2, 0.01, seed)
}

///
/// For posterior sampling debugging
///
pub fn tandem_repeat_small(end_length: usize) -> (Genome, usize) {
    let prefix = generate(end_length, 2);
    let suffix = generate(end_length, 0);

    let unit = b"TAGGACAAGC".to_vec();
    let unit_mut_cg = b"TAGGAGAAGC".to_vec();
    let unit_ins_t = b"TAGGTACAAGC".to_vec();
    let unit_mut_ac = b"TAGGCCAAGC".to_vec();
    let unit_del_a = b"TAGGCAAGC".to_vec();
    let unit_del_g = b"TAGGACAAC".to_vec();

    let hap_0 = StyledSequence::linear(
        [
            prefix.clone(),
            // 50 units,
            tandem_repeat(&unit, 50),
            suffix.clone(),
        ]
        .concat(),
    );
    let hap_1 = StyledSequence::linear(
        [
            prefix,
            // 50 units, 4+20+1+10+10=45 original units and 5 mutated units
            tandem_repeat(&unit, 4),
            unit_mut_cg,
            tandem_repeat(&unit, 20),
            unit_ins_t,
            unit_mut_ac,
            tandem_repeat(&unit, 1),
            unit_del_a,
            tandem_repeat(&unit, 10),
            unit_del_g,
            tandem_repeat(&unit, 10),
            suffix,
        ]
        .concat(),
    );
    let genome_size = hap_0.len() + hap_1.len();

    (vec![hap_0, hap_1], genome_size)
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
    #[test]
    fn genome_tandem_repeat_example_manual_mutation() {
        let (g, gs) = tandem_repeat_diploid_example_mut();
        show_genome(&g, gs);
        assert!(g[0] != g[1]);

        let (g, gs) = tandem_repeat_diploid_example_ins();
        show_genome(&g, gs);
        assert!(g[0] != g[1]);
    }
    #[test]
    fn genome_tandem_repeat_500bp_and_small() {
        let (g_a, gs_a) = tandem_repeat_500bp();
        let (g_b, gs_b) = tandem_repeat_small(50);
        show_genome(&g_a, gs_a);
        show_genome(&g_b, gs_b);
        assert_eq!(g_a[0], g_b[0]);
        for i in 0..(g_b[1].len()) {
            println!(
                "{} {} {}",
                g_a[1].seq()[i],
                g_b[1].seq()[i],
                g_a[1].seq()[i] == g_b[1].seq()[i]
            );
        }
    }
}
