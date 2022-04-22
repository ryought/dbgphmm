#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{sequence_to_string, Genome, Reads, Seq, Sequence};
    use crate::dbg::{Dbg, HashDbg, SimpleDbg};
    use crate::em::compression::{compression, compression_step, compression_with_depths};
    use crate::em::infer;
    use crate::em::scheduler::SchedulerType1;
    use crate::graph::genome_graph::{GenomeGraph, ReadProfile};
    use crate::graph::seq_graph::SeqGraph;
    use crate::hmmv2::params::PHMMParams;
    use crate::hmmv2::sample::{ReadAmount, SampleProfile, StartPoints};
    use crate::kmer::VecKmer;
    use crate::random_seq::{generate, random_mutation, tandem_repeat, MutationProfile};

    fn generate_genome() -> (Genome, usize) {
        let unit_size = 20;
        let n_unit = 20;
        let genome_size = unit_size * n_unit * 2;
        let unit = generate(unit_size, 10);
        let tandem_repeat = tandem_repeat(&unit, n_unit);
        let (hap_a, _) = random_mutation(&tandem_repeat, MutationProfile::uniform(0.1), 0);
        let (hap_b, _) = random_mutation(&hap_a, MutationProfile::uniform(0.01), 2);
        println!("{}", tandem_repeat.to_str());
        println!("{}", hap_a.to_str());
        println!("{}", hap_b.to_str());
        (vec![hap_a, hap_b], genome_size)
    }

    fn generate_dataset() -> (
        Genome,
        Reads,
        SimpleDbg<VecKmer>,
        SimpleDbg<VecKmer>,
        SimpleDbg<VecKmer>,
    ) {
        // (1) generate genome and reads
        println!("generating genome");
        let (genome, genome_size) = generate_genome();
        println!("genome hap_a: {}", sequence_to_string(&genome[0]));
        println!("genome hap_b: {}", sequence_to_string(&genome[1]));

        println!("generating reads");
        let g = GenomeGraph::from_seqs(&genome);
        let profile = ReadProfile {
            has_revcomp: true,
            sample_profile: SampleProfile {
                read_amount: ReadAmount::TotalBases(genome_size * 10),
                seed: 111,
                length: 1000,
                start_points: StartPoints::AllStartPoints,
                endable: false,
            },
            phmm_params: PHMMParams::default(),
        };
        let pos_reads = g.sample_positioned_reads(&profile);
        for pos_read in pos_reads.iter() {
            println!("{}", pos_read);
        }
        g.show_coverage(&pos_reads);
        let reads = pos_reads.to_reads(true);

        let k: usize = 8;
        let dbg_raw: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(k, &reads);

        // (4) compare with true dbg
        let dbg_true_init: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(k, &genome);

        // (5) true k=50 (read length)
        let dbg_true: SimpleDbg<VecKmer> = SimpleDbg::from_seqs(100, &genome);

        (genome, reads, dbg_raw, dbg_true_init, dbg_true)
    }

    #[test]
    fn e2e_tandem_repeat() {
        let (genome, reads, dbg_raw, dbg_true_init, dbg_true) = generate_dataset();

        let scheduler = SchedulerType1::new(8, 100, 10.0);
        let dbg_infer = infer(&dbg_raw, &reads, &PHMMParams::default(), &scheduler, 5);

        println!("{} {}", dbg_true_init.n_traverse_choices(), dbg_true_init);
        println!("{} {}", dbg_true.n_traverse_choices(), dbg_true);
        println!("{} {}", dbg_infer.n_traverse_choices(), dbg_infer);

        let r = dbg_infer.compare(&dbg_true);
        println!("{:?}", r);
    }
}
