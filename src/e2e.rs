//!
//! End-to-end test data generation functions
//!
//! * generate genome
//! * generate reads
//!
use crate::common::{sequence_to_string, Genome, Reads, Seq, Sequence};
use crate::dbg::{Dbg, HashDbg, SimpleDbg};
use crate::graph::genome_graph::{GenomeGraph, ReadProfile};
use crate::hmmv2::params::PHMMParams;
use crate::hmmv2::sample::{ReadAmount, ReadLength, SampleProfile, StartPoints};
use crate::kmer::VecKmer;
use serde::Serialize;
use serde_with::{serde_as, DisplayFromStr};

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

///
/// Dataset collection struct
///
/// * Dataset.reads
/// * Dataset.phmm_params
/// * Dataset.dbg_raw
/// * Dataset.dbg_true_init
/// * Dataset.dbg_true
///
// #[serde_as]
// #[derive(Clone, Serialize)]
pub struct Experiment {
    ///
    /// genome sequence
    ///
    // #[serde_as(as = "Vec<DisplayFromStr>")]
    pub genome: Genome,
    ///
    /// size of the genome
    ///
    pub genome_size: usize,
    ///
    /// sampled reads
    ///
    pub reads: Reads,
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
    ///
    /// show reads
    ///
    pub fn show_reads(&self) {
        for (i, read) in self.reads.iter().enumerate() {
            println!("read#{}\t{}", i, read.to_str());
        }
    }
}

pub fn generate_dataset(
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

    let dbg_raw: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(k_init, &reads);

    // (4) compare with true dbg
    let dbg_true_init: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(k_init, &genome);

    // (5) true k=50 (read length)
    let dbg_true: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(k_target, &genome);

    Experiment {
        genome,
        genome_size,
        reads,
        phmm_params,
        dbg_raw,
        dbg_true_init,
        dbg_true,
    }
}

///
/// will deprecate
///
pub fn generate_full_length_dataset(
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
        genome,
        genome_size,
        reads,
        phmm_params,
        dbg_raw,
        dbg_true_init,
        dbg_true,
    }
}
