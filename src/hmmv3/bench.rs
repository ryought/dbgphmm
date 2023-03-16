//!
//! HMM speed test mock
//!
//! To run test
//!
//! ```text
//! cargo bench -- hmmv2
//! cargo bench -- hmmv2 --nocapture
//! ```
//!

#[cfg(test)]
mod tests {
    use crate::dbg::mocks::mock_random_with_seq;
    use crate::e2e;
    use crate::e2e::Experiment;
    use crate::hmmv2::common::PModel;
    use crate::hmmv2::params::PHMMParams;
    use crate::hmmv2::result::PHMMResultLike;
    use std::time::{Duration, Instant};
    use test::Bencher;

    #[bench]
    fn bench_forward(b: &mut Bencher) {
        let mut param = PHMMParams::default();
        let (dbg, seq) = mock_random_with_seq(8, 100);
        let phmm = dbg.to_phmm(param);
        let r = &seq[20..60];
        b.iter(|| {
            let r1 = phmm.forward(&r);
        });
    }

    // simple genome test
    fn simple_genome() -> (Experiment, PModel) {
        let experiment = e2e::generate_simple_genome_mock();
        let param = experiment.phmm_params;
        let phmm = experiment.dbg_raw.to_phmm(param);
        (experiment, phmm)
    }
    #[bench]
    fn bench_simple_genome_dense_full(b: &mut Bencher) {
        let (experiment, phmm) = simple_genome();
        b.iter(|| {
            let forward = phmm.forward(&experiment.reads()[0].as_ref());
        });
    }
    #[bench]
    fn bench_simple_genome_sparse(b: &mut Bencher) {
        let (experiment, phmm) = simple_genome();
        b.iter(|| {
            let forward = phmm.forward_sparse(&experiment.reads()[0].as_ref());
        });
    }
    #[bench]
    fn bench_simple_genome_dense_with_hint(b: &mut Bencher) {
        let (experiment, phmm) = simple_genome();
        let read = experiment.reads()[0].as_ref();
        let o = phmm.run(read);
        let hint = o.to_hint(10);
        b.iter(|| {
            let forward = phmm.forward_with_hint(read, &hint);
        });
    }
}
