#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::sequence_to_string;
    use crate::dbg::mocks::mock_random_with_seq;
    use crate::hmmv2::params::PHMMParams;
    use crate::hmmv2::result::PHMMResultLike;
    use crate::random_seq::generate;
    #[test]
    fn hmmv2_sparse_and_dense() {
        let mut param = PHMMParams::default();
        // param.n_warmup = 2;
        // param.n_active_nodes = 2;
        let (dbg, seq) = mock_random_with_seq(8, 100);
        let phmm = dbg.to_phmm(param);
        println!("{:?}", sequence_to_string(&seq));

        // let r = generate(100, 1);
        let r = &seq[20..60];
        let r2 = &seq[30..70];
        let r1 = phmm.backward(&r);
        let r2 = phmm.backward_sparse(&r2);
        for i in 0..r1.n_emissions() {
            let d = r1.table(i).diff(&r2.table(i));
            println!("{} {}", i, d);
        }
    }
}
