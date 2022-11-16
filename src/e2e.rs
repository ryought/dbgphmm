//!
//! End-to-end test data generation functions
//!
//! * generate genome
//! * generate reads
//!
use crate::common::{sequence_to_string, Genome, Reads, Seq, Sequence, StyledSequence};
use crate::dbg::{Dbg, HashDbg, SimpleDbg};
use crate::genome;
use crate::graph::genome_graph::{GenomeGraph, ReadProfile};
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::sample::{ReadAmount, ReadLength, SampleProfile, StartPoints};
use crate::kmer::VecKmer;
use serde::{Deserialize, Serialize};

///
/// Read types
///
pub enum ReadType {
    FullLength,
    FixedSizeFragment,
    Fragment,
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
    reads: Reads,
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
    pub fn reads(&self) -> &Reads {
        &self.reads
    }
    pub fn params(&self) -> PHMMParams {
        self.phmm_params
    }
}

///
/// Dataset collection struct
///
/// * Dataset.reads
/// * Dataset.phmm_params
/// * Dataset.dbg_raw
/// * Dataset.dbg_true_init
/// * Dataset.dbg_true
///
#[derive(Clone)]
pub struct Experiment {
    pub dataset: Dataset,
    ///
    /// Profile HMM parameters
    ///
    pub phmm_params: PHMMParams,
    ///
    /// raw read-dbg k=k_init
    ///
    pub dbg_raw: SimpleDbg<VecKmer>,
    ///
    /// true-dbg k=k_init
    ///
    pub dbg_true_init: SimpleDbg<VecKmer>,
    ///
    /// true-dbg k=k_target
    ///
    pub dbg_true: SimpleDbg<VecKmer>,
}

impl Experiment {
    pub fn genome(&self) -> &Genome {
        &self.dataset.genome
    }
    pub fn genome_size(&self) -> usize {
        self.dataset.genome_size
    }
    pub fn reads(&self) -> &Reads {
        &self.dataset.reads
    }
    ///
    /// show reads
    ///
    pub fn show_reads(&self) {
        for (i, read) in self.reads().iter().enumerate() {
            println!("read#{}\t{}", i, read.to_str());
        }
    }
}

pub fn generate_dataset(
    genome: Genome,
    genome_size: usize,
    read_seed: u64,
    coverage: usize,
    read_length: usize,
    read_type: ReadType,
    phmm_params: PHMMParams,
) -> Dataset {
    let g = GenomeGraph::from_styled_seqs(&genome);
    let profile = match read_type {
        ReadType::Fragment => ReadProfile {
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
        ReadType::FullLength => ReadProfile {
            has_revcomp: false,
            sample_profile: SampleProfile {
                read_amount: ReadAmount::Count(coverage),
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
    let reads = pos_reads.to_reads(true);

    Dataset {
        genome,
        genome_size,
        reads,
        phmm_params,
    }
}

pub fn generate_experiment(
    genome: Genome,
    genome_size: usize,
    read_seed: u64,
    phmm_params: PHMMParams,
    coverage: usize,
    read_length: usize,
    read_type: ReadType,
    k_init: usize,
    k_target: usize,
) -> Experiment {
    let dataset = generate_dataset(
        genome,
        genome_size,
        read_seed,
        coverage,
        read_length,
        read_type,
        phmm_params,
    );

    let dbg_raw: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(k_init, dataset.reads());

    // (4) compare with true dbg
    let dbg_true_init: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(k_init, dataset.genome());

    // (5) true k=50 (read length)
    let dbg_true: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(k_target, dataset.genome());

    Experiment {
        dataset,
        phmm_params,
        dbg_raw,
        dbg_true_init,
        dbg_true,
    }
}

///
/// will deprecate
///
pub fn generate_full_length_experiment(
    genome: Genome,
    genome_size: usize,
    read_seed: u64,
    phmm_params: PHMMParams,
    coverage: usize,
) -> Experiment {
    let g = GenomeGraph::from_seqs(&genome);
    let profile = ReadProfile {
        has_revcomp: true,
        sample_profile: SampleProfile {
            read_amount: ReadAmount::TotalBases(genome_size * coverage),
            seed: read_seed,
            length: ReadLength::StateCount(1000),
            start_points: StartPoints::AllStartPoints,
        },
        phmm_params: phmm_params.clone(),
    };
    let pos_reads = g.sample_positioned_reads(&profile);
    g.show_coverage(&pos_reads);
    let reads = pos_reads.to_reads(true);
    for read in reads.iter() {
        println!("{}", read.to_str());
    }

    let k: usize = 8;
    let dbg_raw: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(k, &reads);

    // (4) compare with true dbg
    let dbg_true_init: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(k, &genome);

    // (5) true k=50 (read length)
    let dbg_true: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(100, &genome);

    Experiment {
        dataset: Dataset {
            genome,
            genome_size,
            reads,
            phmm_params,
        },
        phmm_params,
        dbg_raw,
        dbg_true_init,
        dbg_true,
    }
}

///
/// Easy toy example
/// * 100bp simple genome
/// * p=0.1% 20x full length reads
///
pub fn generate_simple_genome_mock() -> Experiment {
    let (genome, genome_size) = genome::simple(100, 5);
    let param = PHMMParams::uniform(0.001);
    generate_experiment(
        genome,
        genome_size,
        0,
        param,
        20, // coverage is 20x
        2000,
        ReadType::FullLength,
        40,
        40,
    )
}

///
/// 1000bp tandem repeat example
///
pub fn generate_small_tandem_repeat() -> Dataset {
    let (genome, genome_size) = genome::tandem_repeat_haploid(100, 10, 0.01, 0, 0);
    let param = PHMMParams::uniform(0.001);
    generate_dataset(
        genome,
        genome_size,
        0,
        10, // coverage is 20x
        2000,
        ReadType::FullLength,
        param,
    )
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn e2e_dataset_serialize_test() {
        let d = Dataset {
            genome: vec![StyledSequence::linear(b"ATCGTTCTTC".to_vec())],
            genome_size: 10,
            reads: Reads::from(vec![b"ATCGT".to_vec(), b"TTTCG".to_vec()]),
            phmm_params: PHMMParams::default(),
        };
        let json = serde_json::to_string(&d).unwrap();
        println!("{}", json);
        assert_eq!(d, serde_json::from_str(&json).unwrap());
    }
}
