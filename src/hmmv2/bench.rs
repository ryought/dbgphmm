//!
//! HMM speed test mock
//!

#[cfg(test)]
mod tests {
    use crate::dbg::mocks::mock_random_with_seq;
    use crate::e2e;
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

    #[bench]
    fn bench_tandem_repeat(b: &mut Bencher) {
        let experiment = e2e::generate_simple_genome_mock();
        let param = experiment.phmm_params;
        let phmm = experiment.dbg_raw.to_phmm(param);

        b.iter(|| {
            let nf = phmm.to_node_freqs(experiment.reads());
        });
    }

    #[test]
    fn bench_tandem_repeat_once() {
        // only measureing the single run
        let experiment = e2e::generate_simple_genome_mock();
        let param = experiment.phmm_params;
        let phmm = experiment.dbg_raw.to_phmm(param);

        let start = Instant::now();
        let nf = phmm.to_node_freqs(experiment.reads());
        let duration = start.elapsed();
        println!("{} ms", duration.as_millis());
    }
}
