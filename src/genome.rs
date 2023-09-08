//!
//! mock genome generation functions
//!
//! * `simple`
//! * `simple_diploid`
//! * `tandem_repeat`
//! * `tandem_repeat_diploid`
//!
use crate::common::{collection::sanitize_bases, SeqStyle, StyledSequence};
use crate::random_seq::{
    generate, join, random_mutation, random_mutation_with_rng, tandem_repeat, MutationProfile,
};
use rand::prelude::*;
use rand_xoshiro::Xoshiro256PlusPlus;
use serde::{Deserialize, Serialize};

/// Genome, the collection of sequences.
///
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Genome(Vec<StyledSequence>);

impl Genome {
    ///
    /// Constructor of Genome
    ///
    pub fn new(seqs: Vec<StyledSequence>) -> Self {
        Genome(seqs)
    }
    /// The number of haplotypes in the genome
    pub fn len(&self) -> usize {
        self.0.len()
    }
    /// The number of linear haplotypes in the genome
    pub fn n_linear_haplotypes(&self) -> usize {
        self.0
            .iter()
            .filter(|h| h.style() == SeqStyle::Linear)
            .count()
    }
    /// Total size (base count) of the genome
    pub fn genome_size(&self) -> usize {
        self.0.iter().map(|seq| seq.len()).sum()
    }
    ///
    /// Print genome as:
    ///
    /// ```text
    /// # genome[0]=ATCGATCGT
    /// # genome[1]=ATCGATCGT
    /// ```
    ///
    pub fn show(&self) {
        for i in 0..self.len() {
            println!("# genome[{}]={}, len={}", i, self[i], self[i].len());
        }
    }
    /// Save genome as FASTA file
    pub fn to_fasta<P: AsRef<std::path::Path>>(&self, path: P) -> std::io::Result<()> {
        let file = std::fs::File::create(path)?;
        self.to_fasta_writer(file)
    }
    /// Save genome as FASTA using `bio::io::fasta::Writer`
    pub fn to_fasta_writer<W: std::io::Write>(&self, writer: W) -> std::io::Result<()> {
        let mut writer = bio::io::fasta::Writer::new(writer);
        for (i, g) in self.into_iter().enumerate() {
            writer.write(&format!("g{}", i), Some(&g.style().to_string()), g.seq())?;
        }
        Ok(())
    }
    /// Parse FASTA file into genome
    pub fn from_fasta<P: AsRef<std::path::Path>>(path: P) -> std::io::Result<Self> {
        let file = std::fs::File::open(path)?;
        Self::from_fasta_reader(file)
    }
    /// Parse FASTA into genome using `bio::io::fasta::Reader`
    ///
    /// If the description string is "C", the haplotype is regarded as circular string
    /// otherwise, it will be linear.
    pub fn from_fasta_reader<R: std::io::Read>(reader: R) -> std::io::Result<Self> {
        let reader = bio::io::fasta::Reader::new(reader);
        let haps = reader
            .records()
            .map(|r| {
                let record = r.unwrap();
                let seq = sanitize_bases(record.seq());
                let style = if record.desc().is_some_and(|desc| desc == "C") {
                    SeqStyle::Circular
                } else {
                    SeqStyle::Linear
                };
                StyledSequence::new(seq, style)
            })
            .collect();
        Ok(Genome(haps))
    }
}

impl std::ops::Index<usize> for Genome {
    type Output = StyledSequence;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl std::ops::IndexMut<usize> for Genome {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<'a> IntoIterator for &'a Genome {
    type Item = &'a StyledSequence;
    type IntoIter = std::slice::Iter<'a, StyledSequence>;
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

///
/// simple random haploid genome
///
/// ## Parameters
///
/// * genome_size: size of the genome
/// * seed: seed of the sequence generation
///
pub fn simple(genome_size: usize, seed: u64) -> Genome {
    Genome::new(vec![StyledSequence::linear(generate(genome_size, seed))])
}

///
/// simple random haploid circular genome
///
pub fn simple_circular(genome_size: usize, seed: u64) -> Genome {
    Genome::new(vec![StyledSequence::circular(generate(genome_size, seed))])
}

///
/// simple diploid genome, with two SNVs manually added.
///
/// will be deprecated
///
pub fn simple_diploid() -> Genome {
    let genome_size_hap = 100;
    let haplotype1 = generate(100, 0);
    let mut haplotype2 = haplotype1.clone();
    haplotype2[30] = b'C';
    haplotype2[80] = b'T';
    Genome::new(vec![
        StyledSequence::linear(haplotype1),
        StyledSequence::linear(haplotype2),
    ])
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
pub fn diploid(hap_size: usize, hap_seed: u64, div_rate: f64, div_seed: u64) -> Genome {
    let hap_a = generate(hap_size, hap_seed);
    let (hap_b, ops) = random_mutation(&hap_a, MutationProfile::uniform(div_rate), div_seed);
    // println!("ops={:?}", ops);
    Genome::new(vec![
        StyledSequence::linear(hap_a),
        StyledSequence::linear(hap_b),
    ])
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
) -> Genome {
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
) -> Genome {
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
    Genome::new(vec![StyledSequence::linear(hap_a)])
}

pub fn tandem_repeat_diploid(
    unit_size: usize,
    n_unit: usize,
    divergence_init: f64,
    unit_seed: u64,
    hap_seed: u64,
    divergence_between_haplotypes: f64,
    div_seed: u64,
) -> Genome {
    let mut hap = tandem_repeat_haploid(unit_size, n_unit, divergence_init, unit_seed, hap_seed);
    let hap_a = hap[0].clone();
    let (hap_b_seq, ops) = random_mutation(
        &hap_a.seq(),
        MutationProfile::uniform(divergence_between_haplotypes),
        div_seed,
    );
    let hap_b = StyledSequence::linear(hap_b_seq);
    Genome::new(vec![hap_a, hap_b])
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
) -> Genome {
    let mut hap = tandem_repeat_haploid_with_unique_ends(
        unit_size,
        n_unit,
        divergence_init,
        unit_seed,
        hap_seed,
        end_length,
    );
    let mut seqs = Vec::new();
    let hap_a = hap[0].clone();
    seqs.push(hap_a.clone());

    let mut rng = Xoshiro256PlusPlus::seed_from_u64(div_seed);
    for _ in 1..n_haplotypes {
        let (hap_b_seq, ops) = random_mutation_with_rng(
            &hap_a.seq(),
            MutationProfile::uniform(divergence_between_haplotypes),
            &mut rng,
        );
        let hap_b = StyledSequence::linear(hap_b_seq);
        seqs.push(hap_b);
    }
    Genome::new(seqs)
}

///
/// Tandem repeat (polyploid) genome with unique (homozygous) prefix/suffix.
///
pub fn tandem_repeat_polyploid_with_unique_homo_ends(
    unit_size: usize,
    n_unit: usize,
    unit_seed: u64,
    divergence_init: f64,
    div_init_seed: u64,
    end_length: usize,
    n_haplotypes: usize,
    divergence_between_haplotypes: f64,
    div_seed: u64,
) -> Genome {
    // base tandem repeats
    let unit = generate(unit_size, unit_seed);
    let tandem_repeat = tandem_repeat(&unit, n_unit);

    // initial divergence
    let (tandem_repeat, ops) = random_mutation(
        &tandem_repeat,
        MutationProfile::uniform(divergence_init),
        div_init_seed,
    );
    // println!("[genome] ops hap[0] {:?}", ops);

    // unique ends
    let prefix = generate(end_length, unit_seed.wrapping_add(1));
    let suffix = generate(end_length, unit_seed.wrapping_sub(1));

    let mut genome = Vec::new();

    // add first (unmutated) haplotype
    let hap = join(prefix.clone(), join(tandem_repeat.clone(), suffix.clone()));
    genome.push(StyledSequence::linear(hap));

    // add other mutated haplotypes
    let mut rng = Xoshiro256PlusPlus::seed_from_u64(div_seed);
    for i in 1..n_haplotypes {
        let (tandem_repeat_mut, ops) = random_mutation_with_rng(
            &tandem_repeat,
            MutationProfile::uniform(divergence_between_haplotypes),
            &mut rng,
        );
        println!("[genome] ops hap[{}] {:?}", i, ops);
        let hap = join(prefix.clone(), join(tandem_repeat_mut, suffix.clone()));
        genome.push(StyledSequence::linear(hap));
    }

    Genome::new(genome)
}

pub fn tandem_repeat_diploid_example_mut() -> Genome {
    let mut genome = tandem_repeat_polyploid_with_unique_ends(50, 20, 0.00, 0, 0, 50, 2, 0.00, 0);
    genome[1].seq[75] = b'A';
    genome
}

pub fn tandem_repeat_diploid_example_ins() -> Genome {
    let mut genome = tandem_repeat_polyploid_with_unique_ends(50, 20, 0.00, 0, 0, 50, 2, 0.00, 0);
    genome[1].seq.insert(75, b'C');
    genome
}

pub fn tandem_repeat_500bp() -> Genome {
    let seed = 1;
    tandem_repeat_polyploid_with_unique_ends(10, 50, 0.0, seed, seed, 50, 2, 0.01, seed)
}

///
/// For posterior sampling debugging
///
/// original:
/// let (a, b, c, d) = (4, 20, 10, 10);
/// test:
/// let (a, b, c, d) = (1, 2, 1, 20);
///
pub fn tandem_repeat_small(
    end_length: usize,
    a: usize,
    b: usize,
    c: usize,
    d: usize,
    mut_cg: bool,
    ins_t: bool,
    mut_ac: bool,
    del_a: bool,
    del_g: bool,
) -> Genome {
    let prefix = generate(end_length, 2);
    let suffix = generate(end_length, 0);

    let unit = b"TAGGACAAGC".to_vec();
    let unit_mut_cg = b"TAGGAGAAGC".to_vec();
    let unit_ins_t = b"TAGGTACAAGC".to_vec();
    let unit_mut_ac = b"TAGGCCAAGC".to_vec();
    let unit_del_a = b"TAGGCAAGC".to_vec(); // missing
    let unit_del_g = b"TAGGACAAC".to_vec();

    // original
    let n = a + b + c + d + 6;
    let hap_0 = StyledSequence::linear(
        [
            prefix.clone(),
            // 50 units,
            tandem_repeat(&unit, n),
            suffix.clone(),
        ]
        .concat(),
    );
    let hap_1 = StyledSequence::linear(
        [
            prefix,
            // 50 units, 4+20+1+10+10=45 original units and 5 mutated units
            tandem_repeat(&unit, a),
            if mut_cg { unit_mut_cg } else { unit.clone() },
            tandem_repeat(&unit, b),
            if ins_t { unit_ins_t } else { unit.clone() },
            if mut_ac { unit_mut_ac } else { unit.clone() },
            tandem_repeat(&unit, 1),
            if del_a { unit_del_a } else { unit.clone() }, // missing
            tandem_repeat(&unit, c),
            if del_g { unit_del_g } else { unit.clone() }, // missing
            tandem_repeat(&unit, d),
            suffix,
        ]
        .concat(),
    );

    Genome::new(vec![hap_0, hap_1])
}

/// Test case genome `n`
/// 500bp unit `n`-times repeat, 300bp unique ends
/// diploid with H0=2%, H=2%
pub fn n(n: usize) -> Genome {
    tandem_repeat_polyploid_with_unique_homo_ends(500, n, 0, 0.02, 1, 300, 2, 0.02, 0)
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn e2e_genome_generation() {
        let g = simple(100, 0);
        g.show();
        assert_eq!(g.genome_size(), 100);
        assert_eq!(g.len(), 1);
        assert_eq!(g, Genome::new(vec![StyledSequence::linear(b"CCAATTCACAAAAACCACACCTTGGCCAAGGTATCGTATCTTGTTGTTGTATGTGAAAGGGGCCCTAAGATCTGTAGCCACCATGGCTAGGGTCAAATCT".to_vec())]));

        let g = simple(50, 5);
        g.show();
        assert_eq!(g.genome_size(), 50);
        assert_eq!(g.len(), 1);
        assert_eq!(
            g,
            Genome::new(vec![StyledSequence::linear(
                b"CGAAGATGAGAAACCCGGGAGTCGATATATTCAAACAAACGGGCGCTCCT".to_vec()
            )])
        );

        let g = simple_diploid();
        g.show();
        assert_eq!(g.genome_size(), 200);
        assert_eq!(g.len(), 2);
        assert_eq!(
            g,
            Genome::new(vec![
                StyledSequence::linear(b"CCAATTCACAAAAACCACACCTTGGCCAAGGTATCGTATCTTGTTGTTGTATGTGAAAGGGGCCCTAAGATCTGTAGCCACCATGGCTAGGGTCAAATCT".to_vec()),
                StyledSequence::linear(b"CCAATTCACAAAAACCACACCTTGGCCAAGCTATCGTATCTTGTTGTTGTATGTGAAAGGGGCCCTAAGATCTGTAGCCATCATGGCTAGGGTCAAATCT".to_vec()),
            ])
        );

        let g = diploid(50, 5, 0.1, 0);
        g.show();
        assert_eq!(g.genome_size(), 101);
        assert_eq!(g.len(), 2);
        assert_eq!(
            g,
            Genome::new(vec![
                StyledSequence::linear(
                    b"CGAAGATGAGAAACCCGGGAGTCGATATATTCAAACAAACGGGCGCTCCT".to_vec()
                ),
                StyledSequence::linear(
                    b"CGAAGAATGAGAAACACCGGAGTCGTATATTCCAAACAAACGGGCGCTCCT".to_vec()
                )
            ])
        );

        let g = tandem_repeat_haploid(10, 5, 0.0, 0, 0);
        g.show();
        assert_eq!(g.genome_size(), 50);
        assert_eq!(g.len(), 1);
        assert_eq!(
            g,
            Genome::new(vec![StyledSequence::linear(
                b"CCAATTCACACCAATTCACACCAATTCACACCAATTCACACCAATTCACA".to_vec()
            )])
        );

        let g = tandem_repeat_haploid(10, 5, 0.1, 55, 33);
        g.show();
        assert_eq!(g.genome_size(), 52);
        assert_eq!(g.len(), 1);
        assert_eq!(
            g,
            Genome::new(vec![StyledSequence::linear(
                b"GTAAATGCGGGTAAATTGCGGGTAAATGCGGCGTAAATGCGGGGAAATCGGG".to_vec()
            )])
        );

        let g = tandem_repeat_diploid(10, 5, 0.1, 55, 33, 0.05, 22);
        g.show();
        assert_eq!(g.genome_size(), 106);
        assert_eq!(g.len(), 2);
        assert_eq!(
            g,
            Genome::new(vec![
                StyledSequence::linear(
                    b"GTAAATGCGGGTAAATTGCGGGTAAATGCGGCGTAAATGCGGGGAAATCGGG".to_vec()
                ),
                StyledSequence::linear(
                    b"GTAAATGCGGGTAACTCTGCGGGTAAATGCGGCGTAAATGCGGGGACAATCGGG".to_vec()
                ),
            ])
        );

        let g = tandem_repeat_haploid_with_unique_ends(10, 5, 0.0, 0, 0, 10);
        g.show();
        assert_eq!(g.genome_size(), 70);
        assert_eq!(g.len(), 1);
        assert_eq!(
            g,
            Genome::new(vec![StyledSequence::linear(
                b"TAGGACAAGCCCAATTCACACCAATTCACACCAATTCACACCAATTCACACCAATTCACACCTCACCTCA".to_vec()
            )])
        );
    }
    #[test]
    fn genome_tandem_repeat_polyploid() {
        // ploidy=1
        let g = tandem_repeat_polyploid_with_unique_ends(10, 5, 0.0, 0, 0, 0, 1, 0.01, 0);
        g.show();
        assert_eq!(
            g,
            Genome::new(vec![StyledSequence::linear(
                b"CCAATTCACACCAATTCACACCAATTCACACCAATTCACACCAATTCACA".to_vec()
            )])
        );

        // ploidy=4
        let g = tandem_repeat_polyploid_with_unique_ends(10, 5, 0.0, 0, 0, 0, 4, 0.01, 0);
        g.show();
        assert_eq!(
            g,
            Genome::new(vec![
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
            ])
        );
    }
    #[test]
    fn genome_tandem_repeat_example_manual_mutation() {
        let g = tandem_repeat_diploid_example_mut();
        g.show();
        assert!(g[0] != g[1]);

        let g = tandem_repeat_diploid_example_ins();
        g.show();
        assert!(g[0] != g[1]);
    }
    #[test]
    fn genome_tandem_repeat_500bp_and_small() {
        let g_a = tandem_repeat_500bp();
        let g_b = tandem_repeat_small(50, 4, 20, 10, 10, true, true, true, true, true);
        g_a.show();
        g_b.show();
        assert_eq!(g_a[0], g_b[0]);
        for i in 0..(g_b[1].len()) {
            println!(
                "{} {} {}",
                g_a[1].seq()[i],
                g_b[1].seq()[i],
                g_a[1].seq()[i] == g_b[1].seq()[i]
            );
        }
        assert_eq!(g_a[1].seq()[..550], g_b[1].seq()[..550]);
    }
    #[test]
    fn genome_tandem_repeat_unique_homo_ends() {
        let g = tandem_repeat_polyploid_with_unique_homo_ends(10, 5, 0, 0.0, 0, 10, 4, 0.05, 0);
        g.show();
        assert_eq!(
            g,
            Genome::new(vec![
                StyledSequence::linear(
                    b"TAGGACAAGCCCAATTCACACCAATTCACACCAATTCACACCAATTCACACCAATTCACACCTCACCTCA"
                        .to_vec()
                ),
                StyledSequence::linear(
                    b"TAGGACAAGCCCAATTCACACCAAATGACACCAATCACACCAATTCACACCAATTCACACCTCACCTCA"
                        .to_vec()
                ),
                StyledSequence::linear(
                    b"TAGGACAAGCCCAATATCCACACCAATTCACACCAATTCACACCCAATTCACACCAATTCACACCTCACCTCA"
                        .to_vec()
                ),
                StyledSequence::linear(
                    b"TAGGACAAGCCCAATTCACAGCAATTCACCCAATTCACACCTATTCACACCAATTCACACCTCACCTCA"
                        .to_vec()
                ),
            ])
        );
    }
    #[test]
    fn genome_fasta() {
        let genome = Genome::new(vec![
            StyledSequence::linear(b"ATCGGCT".to_vec()),
            StyledSequence::linear(b"TTCGCC".to_vec()),
            StyledSequence::circular(b"AAATAA".to_vec()),
        ]);

        // save
        let mut file: Vec<u8> = Vec::new();
        genome.to_fasta_writer(&mut file);
        let s = String::from_utf8(file.clone()).unwrap();
        println!("{}", s);
        assert_eq!(s, ">g0 L\nATCGGCT\n>g1 L\nTTCGCC\n>g2 C\nAAATAA\n");

        // load
        let genome2 = Genome::from_fasta_reader(file.as_slice()).unwrap();
        genome2.show();
        assert_eq!(genome2, genome);
    }
}
