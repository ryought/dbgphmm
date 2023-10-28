//!
//! End-to-end test data generation functions
//!
//! * generate genome
//! * generate reads
//!
use crate::common::{Genome, PositionedReads, Seq};
use crate::genome;
use crate::graph::genome_graph::{GenomeGraph, ReadProfile};
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::sample::{ReadAmount, ReadLength, SampleProfile, StartPoints};
use crate::utils::spaces;
use serde::{Deserialize, Serialize};
use std::io::Write;

///
/// Read types
///
pub enum ReadType {
    FullLengthForHaploid,
    FullLength,
    FixedSizeFragment,
    FragmentWithRevComp,
    FullLengthWithRevComp,
}

///
/// Dataset configuration
///
pub struct DatasetConfig {}

#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct Dataset {
    ///
    /// genome sequence
    ///
    genome: Genome,
    ///
    /// size of the genome
    ///
    genome_size: usize,
    ///
    /// sampled reads
    ///
    reads: PositionedReads,
    ///
    /// Profile HMM parameters used in read sampling
    ///
    phmm_params: PHMMParams,
    // ///
    // /// dataset label
    // ///
    // pub label: String,
}

impl Dataset {
    pub fn genome(&self) -> &Genome {
        &self.genome
    }
    pub fn genome_size(&self) -> usize {
        self.genome_size
    }
    pub fn reads(&self) -> &PositionedReads {
        &self.reads
    }
    pub fn params(&self) -> PHMMParams {
        self.phmm_params
    }
    ///
    /// estimate coverage of reads by (total_bases_in_reads) / (true_genome_size)
    ///
    pub fn coverage(&self) -> f64 {
        self.reads().coverage(self.genome_size())
    }
    ///
    /// Average of read lengths
    ///
    pub fn average_read_length(&self) -> usize {
        self.reads().average_length()
    }
    ///
    /// show reads
    ///
    pub fn show_reads(&self) {
        self.reads().show_reads()
    }
    pub fn show_genome(&self) {
        self.genome().show()
    }
    ///
    ///
    ///
    pub fn show_reads_with_genome(&self) {
        let _n = self.reads().len();
        let n_hap = self.genome().len();
        let mut pos: Vec<Vec<usize>> = vec![vec![]; n_hap];
        for read_id in 0..(self.reads().len()) {
            // haplotype id that generated the read
            let hap_id = self.reads()[read_id].origin_node().index();
            pos[hap_id].push(read_id);
        }
        for hap_id in 0..n_hap {
            pos[hap_id].sort_by_key(|&read_id| self.reads()[read_id].origin_pos())
        }

        for hap_id in 0..n_hap {
            for (i, &read_id) in pos[hap_id].iter().enumerate() {
                if i % 10 == 0 {
                    println!("# g[{}]\n# {}", hap_id, self.genome()[hap_id]);
                }
                let read = &self.reads()[read_id];
                let offset = 2 + read.origin_pos();
                let spaces = spaces(offset);
                // let aligned = read.to_aligned_str();
                // println!(
                //     "{}r[{}]\n{}{}\n{}{}",
                //     spaces, read_id, spaces, aligned[0], spaces, aligned[1],
                // );
                println!("# {}{} r[{}]", spaces, read.seq().to_str(), read_id);
            }
        }
    }
    pub fn to_json_file<P: AsRef<std::path::Path>>(&self, path: P) {
        let mut file = std::fs::File::create(path).unwrap();
        serde_json::to_writer(&mut file, self).unwrap()
    }
    pub fn from_json_file<P: AsRef<std::path::Path>>(path: P) -> Self {
        let mut file = std::fs::File::open(path).unwrap();
        serde_json::from_reader(&mut file).unwrap()
    }
    ///
    /// Dump genome as fasta
    ///
    pub fn to_genome_fasta<P: AsRef<std::path::Path>>(&self, path: P) -> std::io::Result<()> {
        self.genome().to_fasta(path)
    }
    ///
    /// Dump reads as fasta.
    /// Sample position and hint information will be removed.
    ///
    pub fn to_reads_fasta<P: AsRef<std::path::Path>>(&self, path: P) -> std::io::Result<()> {
        self.reads().to_fasta(path)
    }
    ///
    /// Dump SAM file of reads
    ///
    pub fn to_reads_sam<P: AsRef<std::path::Path>>(&self, path: P) -> std::io::Result<()> {
        let mut file = std::fs::File::create(path).unwrap();

        // header
        for (i, g) in self.genome().into_iter().enumerate() {
            writeln!(file, "@SQ\tSN:g{}\tLN:{}", i, g.len())?;
        }

        for (i, r) in self.reads().iter().enumerate() {
            writeln!(file, "{}", r.to_sam_string(&format!("r{}", i)))?;
        }

        Ok(())
    }
}

pub fn generate_dataset(
    genome: Genome,
    read_seed: u64,
    coverage: usize,
    read_length: usize,
    read_type: ReadType,
    phmm_params: PHMMParams,
) -> Dataset {
    let g = GenomeGraph::from_styled_seqs(&genome);
    let genome_size = genome.genome_size();
    let profile = match read_type {
        ReadType::FragmentWithRevComp => ReadProfile {
            has_revcomp: true,
            sample_profile: SampleProfile {
                read_amount: ReadAmount::TotalBases(genome_size * coverage),
                seed: read_seed,
                length: ReadLength::StateCount(read_length),
                start_points: StartPoints::Random,
            },
            phmm_params: phmm_params.clone(),
        },
        ReadType::FixedSizeFragment => ReadProfile {
            has_revcomp: false,
            sample_profile: SampleProfile {
                read_amount: ReadAmount::TotalBases(genome_size * coverage),
                seed: read_seed,
                length: ReadLength::EmitCount(read_length),
                start_points: StartPoints::Random,
            },
            phmm_params: phmm_params.clone(),
        },
        ReadType::FullLengthForHaploid => ReadProfile {
            has_revcomp: false,
            sample_profile: SampleProfile {
                read_amount: ReadAmount::Count(coverage),
                seed: read_seed,
                length: ReadLength::StateCount(read_length),
                start_points: StartPoints::AllStartPoints,
            },
            phmm_params: phmm_params.clone(),
        },
        ReadType::FullLength => ReadProfile {
            has_revcomp: false,
            sample_profile: SampleProfile {
                read_amount: ReadAmount::TotalBases(genome_size * coverage),
                seed: read_seed,
                length: ReadLength::StateCount(read_length),
                start_points: StartPoints::AllStartPoints,
            },
            phmm_params: phmm_params.clone(),
        },
        ReadType::FullLengthWithRevComp => ReadProfile {
            has_revcomp: true,
            sample_profile: SampleProfile {
                read_amount: ReadAmount::TotalBases(genome_size * coverage),
                seed: read_seed,
                length: ReadLength::StateCount(read_length),
                start_points: StartPoints::AllStartPoints,
            },
            phmm_params: phmm_params.clone(),
        },
    };
    let pos_reads = g.sample_positioned_reads(&profile);
    // for read in pos_reads.iter() {
    //     println!("{}", read);
    // }
    // g.show_coverage(&pos_reads);
    // TODO
    // use strand justified read only currently
    let reads = pos_reads.justify_strand();

    Dataset {
        genome,
        genome_size,
        reads,
        phmm_params,
    }
}

///
/// Easy toy example
/// * 200bp simple genome
/// * p=0.1% 20x fragment reads 50bp fixed length reads
///
pub fn generate_simple_genome_fragment_dataset() -> Dataset {
    let genome = genome::simple(200, 5);
    let param = PHMMParams::uniform(0.001);
    generate_dataset(
        genome,
        0,
        20, // coverage (20x)
        50, // length (50bp)
        ReadType::FragmentWithRevComp,
        param,
    )
}

///
/// Easy toy example
/// * 200bp tandem repeat (20bp-unique-prefix + 40bp x 4 + 20bp-unique-suffix) genome
/// * p=0.1% 20x fragment reads 50bp fixed length reads
///
pub fn generate_tandem_repeat_fragment_dataset() -> Dataset {
    let genome = genome::tandem_repeat_haploid_with_unique_ends(40, 4, 0.01, 0, 0, 20);
    let param = PHMMParams::uniform(0.001);
    generate_dataset(
        genome,
        0,
        20, // coverage (20x)
        50, // length (50bp)
        ReadType::FragmentWithRevComp,
        param,
    )
}

///
/// "-c 20 -l 100 -p 0.01 --k-init 12 --k-final 100 -U 50 -N 20 -E 50 -P 2 -D 0.05 -H 0.05 --sigma 100 -m 10 --use-fragment-read";;
///
pub fn generate_difficult_diploid_tandem_repeat_dataset() -> Dataset {
    let genome =
        genome::tandem_repeat_polyploid_with_unique_ends(50, 20, 0.05, 0, 0, 50, 2, 0.05, 0);
    let param = PHMMParams::uniform(0.01);
    generate_dataset(
        genome,
        0,  // read seed
        20, // coverage
        100,
        ReadType::FragmentWithRevComp,
        param,
    )
}

///
///
///
pub fn generate_difficult_diploid_tandem_repeat_dataset_full_length() -> Dataset {
    let genome =
        genome::tandem_repeat_polyploid_with_unique_ends(50, 20, 0.05, 0, 0, 50, 2, 0.05, 0);
    let genome_size = genome.genome_size();
    let param = PHMMParams::uniform(0.01);
    generate_dataset(
        genome,
        0,  // read seed
        20, // coverage
        genome_size * 2,
        ReadType::FullLength,
        param,
    )
}

///
/// Difficult 500bp tandem repeat test case
///
/// U10N50H001S1
///
/// -c 20 -l 100 -p 0.01
/// --k-init 12 --k-final 15
/// -U 10 -N 50 -E 50 -P 2 -D 0.0 -H 0.01 --sigma 100 -d 10 -m 50 -s 1
/// --use-true-end-nodes --start-from-true --use-true-dbg
///
pub fn generate_500bp_case_dataset() -> Dataset {
    let seed = 1;
    let genome = genome::tandem_repeat_500bp();
    let genome_size = genome.genome_size();
    let param = PHMMParams::uniform(0.01);
    generate_dataset(
        genome,
        seed,
        20, // coverage
        genome_size * 2,
        ReadType::FullLength,
        param,
    )
}

///
/// small tandem repeat test case
///
pub fn generate_small_case_dataset(
    a: usize,
    b: usize,
    c: usize,
    d: usize,
    mut_cg: bool,
    ins_t: bool,
    mut_ac: bool,
    del_a: bool,
    del_g: bool,
    p: f64,
) -> Dataset {
    let genome = genome::tandem_repeat_small(20, a, b, c, d, mut_cg, ins_t, mut_ac, del_a, del_g);
    let genome_size = genome.genome_size();
    let seed = 1;
    let param = PHMMParams::uniform(p);
    generate_dataset(
        genome,
        seed,
        20, // coverage
        genome_size * 2,
        ReadType::FullLength,
        param,
    )
}

pub fn generate_tandem_repeat_1kbp() -> Dataset {
    let seed = 1;
    let genome = genome::tandem_repeat_polyploid_with_unique_ends(
        50, 20, 0.0, seed, seed, 50, 2, 0.01, seed,
    );
    let seed = 1;
    let param = PHMMParams::uniform(0.01);
    generate_dataset(
        genome,
        seed,
        20, // coverage
        200,
        ReadType::FixedSizeFragment,
        param,
    )
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ni;
    use crate::common::{PositionedSequence, StyledSequence};
    use crate::graph::genome_graph::GenomeGraphPos;

    #[test]
    fn e2e_dataset_serialize_test() {
        let d = Dataset {
            genome: Genome::new(vec![StyledSequence::linear(b"ATCGTTCTTC".to_vec())]),
            genome_size: 10,
            reads: PositionedReads::from(vec![
                PositionedSequence::new(
                    b"ATCGT".to_vec(),
                    vec![
                        GenomeGraphPos::new_match(ni(0), 0),
                        GenomeGraphPos::new_match(ni(0), 1),
                        GenomeGraphPos::new_match(ni(0), 2),
                        GenomeGraphPos::new_match(ni(0), 3),
                        GenomeGraphPos::new_match(ni(0), 3),
                    ],
                    false,
                ),
                PositionedSequence::new(
                    b"TTTCG".to_vec(),
                    vec![
                        GenomeGraphPos::new_match(ni(1), 10),
                        GenomeGraphPos::new_match(ni(1), 9),
                        GenomeGraphPos::new_match(ni(1), 8),
                        GenomeGraphPos::new_match(ni(1), 6),
                        GenomeGraphPos::new_match(ni(1), 5),
                    ],
                    true,
                ),
            ]),
            phmm_params: PHMMParams::default(),
        };

        // json
        let json = serde_json::to_string(&d).unwrap();
        println!("{:?}", json);
        assert_eq!(d, serde_json::from_str(&json).unwrap());

        // fasta
        d.to_genome_fasta("d.genome.fasta");
        d.to_reads_fasta("d.reads.fasta");
    }
    #[test]
    fn e2e_dataset_read_with_genome_visualization_test() {
        println!("generating genome");
        let genome =
            genome::tandem_repeat_polyploid_with_unique_ends(50, 4, 0.05, 0, 0, 50, 2, 0.05, 0);
        println!("generating reads");
        let param = PHMMParams::uniform(0.01);
        let dataset = generate_dataset(
            genome,
            0,  // read seed
            20, // coverage
            100,
            ReadType::FragmentWithRevComp,
            param,
        );
        dataset.show_reads_with_genome();
    }
}
