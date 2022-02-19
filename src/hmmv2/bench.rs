//!
//! HMM speed test mock
//!

#[cfg(test)]
mod tests {
    use crate::dbg::mocks::mock_random_with_seq;
    use crate::hmmv2::params::PHMMParams;
    use crate::hmmv2::result::PHMMResultLike;
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
}
