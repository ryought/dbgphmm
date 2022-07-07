//!
//! Fragmented read test
//!
//! ## (Ignored) tests
//!
//! * e2e_fragment_full
//!
#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{sequence_to_string, Genome, Reads, Seq, Sequence};
    use crate::dbg::{Dbg, HashDbg, SimpleDbg};
    use crate::e2e::{generate_dataset, Dataset, ReadType};
    use crate::em::compression::{compression, compression_step, compression_with_depths};
    use crate::em::e2e::runner::{benchmark, show_logs};
    use crate::em::TaskLog;
    use crate::genome::simple;
    use crate::graph::genome_graph::{GenomeGraph, ReadProfile};
    use crate::graph::seq_graph::SeqGraph;
    use crate::hmmv2::params::PHMMParams;
    use crate::hmmv2::sample::{ReadAmount, SampleProfile, StartPoints};
    use crate::kmer::VecKmer;
    use crate::random_seq::generate;

    fn generate_e2e_fragment_mock() -> Dataset {
        let (genome, genome_size) = simple(200, 0);
        println!("genome: {}", sequence_to_string(&genome[0]));

        let g = GenomeGraph::from_seqs(&genome);
        let coverage = 10;
        let read_length = 50;
        generate_dataset(
            genome,                //
            genome_size,           //
            11,                    // read_seed
            PHMMParams::default(), //
            coverage,              //
            read_length,           //
            ReadType::Fragment,    //
            8,                     // k_init
            read_length,           // k_target
        )
    }

    #[test]
    fn e2e_fragment_full_benchmark_em_steps() {
        let dataset = generate_e2e_fragment_mock();
        /*
        benchmark_em_steps(
            &dbg_raw,
            &dbg_true,
            &reads,
            &genome,
            &PHMMParams::default(),
            10.0,
        );
        */
        let (dbg, logs) = compression(
            &dataset.dbg_raw,
            &dataset.reads,
            &PHMMParams::default(),
            1.0,
            10,
        );
        show_logs(&[TaskLog::Compression(logs)], &dataset);
    }

    #[test]
    fn e2e_fragment_full() {
        let dataset = generate_e2e_fragment_mock();
        let (dbg_infer, r, _) = benchmark(&dataset, 10.0);

        /*
        let (dbg, _) = compression(&dbg_raw, &reads, &PHMMParams::default(), 1.0, 10);
        let (dbg, logs) = compression_with_depths(
            &dbg_raw,
            &reads,
            &PHMMParams::default(),
            &[
                1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 5.0, 5.0, 8.0, 8.0, 10.0, 10.0,
            ],
        );
        println!("{}", dbg);
        println!("{:?}", logs);
        println!("{}", dbg_true);

        let r = dbg.compare(&dbg_true);
        println!("{:?}", r);
        */
    }
}
