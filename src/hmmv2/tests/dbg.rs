//!
//! Accuracy/speed of MultiDbg PHMM
//!
//! ## Comparison
//! * with/without mapping
//! * non-zero/normal/uniform phmm
//! * adaptive-sparse/normal-sparse/dense
//! * [`test_read_dbg_n_warmup`]
//!     n_warmup fixed/adaptive
//!
//!
//! ## DBG
//! * from genome (without sequencing error)
//! * (TODO) from reads (with error)
//!
//! ## Dataset
//! * unique
//! * small tandem repeat
//! * large tandem repeat
//!

#[cfg(test)]
mod tests {
    //
    // speed benchmarks
    //
    use test::Bencher;

    //
    // common settings
    //
    use crate::e2e::{generate_dataset, Dataset, ReadType};
    use crate::genome;
    use crate::hmmv2::common::PModel;
    use crate::hmmv2::sample::State;
    use crate::hmmv2::table::PHMMTable;
    use crate::multi_dbg::MultiDbg;
    use crate::prelude::*;
    use crate::prob::{lp, p};
    use crate::utils::timer;
    use itertools::Itertools;
    use std::io::prelude::*;

    ///
    /// * Dataset
    /// * initial MultiDbg
    ///
    fn dataset(genome: Genome, k: usize) -> (Dataset, MultiDbg) {
        let dataset = generate_dataset(
            genome,
            0,
            20,                            // 20x coverage
            1000,                          // 1000bp reads
            ReadType::FragmentWithRevComp, // fragment read
            PHMMParams::uniform(0.001),    // 0.1% HiFi error
        );
        let (dbg, t) = timer(|| MultiDbg::create_draft_from_dataset(k, &dataset));
        println!("t={}", t);
        (dataset, dbg)
    }

    /// unit 1kbp x 1 times = 1kB diploid (unique sequence)
    fn u1k() -> Genome {
        genome::tandem_repeat_polyploid_with_unique_homo_ends(1000, 1, 0, 0.0, 0, 50, 2, 0.01, 0)
    }
    /// unit 20bp x 50 times = 1kB diploid
    fn u20() -> Genome {
        genome::tandem_repeat_polyploid_with_unique_homo_ends(20, 50, 0, 0.0, 0, 50, 2, 0.01, 0)
    }
    /// unit 20bp x 200 times = 4kB diploid
    fn u20n200() -> Genome {
        genome::tandem_repeat_polyploid_with_unique_homo_ends(20, 200, 0, 0.02, 0, 300, 2, 0.02, 0)
    }
    /// unit 100bp x 10 times = 1kB diploid
    fn u100() -> Genome {
        genome::tandem_repeat_polyploid_with_unique_homo_ends(100, 10, 0, 0.0, 0, 50, 2, 0.01, 0)
    }

    ///
    /// compare with/without mapping of read dbg
    ///
    fn test_read_dbg_mapping(genome: Genome, k: usize) {
        // create read DBG
        let (dataset, dbg) = dataset(genome, k);
        // dbg.to_gfa_file("read_dbg.gfa");
        let (mapping, t) = timer(|| dbg.generate_mappings(dataset.params(), dataset.reads(), None));
        println!("mapping created in t={}", t);
        let phmm = dbg.to_phmm(dataset.params());

        //
        // (1) full-prob is same with/with-out mapping
        // likelihood of reads
        //
        // 0: without mapping
        let (p0, t0) = timer(|| dbg.to_likelihood(dataset.params(), dataset.reads(), None));
        // 1: with mapping
        let (p1, t1) =
            timer(|| dbg.to_likelihood(dataset.params(), dataset.reads(), Some(&mapping)));
        // 2: without mapping and do not use score only calculation
        let (p2, t2) = timer(|| phmm.to_full_prob_sparse(dataset.reads(), true));
        let diff0 = p0.log_diff(p2);
        let diff1 = p1.log_diff(p2);
        println!(
            "p0={} p1={} p2={} diff0={} diff1={} t0={} t1={} t2={}",
            p0, p1, p2, diff0, diff1, t0, t1, t2
        );
        assert!(diff0 < 0.0001);
        assert!(diff1 < 0.0001);

        // likelihood of genomes TODO
    }
    ///
    /// compare with and without max_ratio
    ///
    fn test_read_dbg_max_ratio(genome: Genome, k: usize) {
        // Create read DBG
        let (dataset, dbg) = dataset(genome, k);

        // Create three patterns of PHMM
        // 0
        let phmm0 = dbg.to_phmm(dataset.params());

        // 1
        let mut param1 = dataset.params();
        param1.n_active_nodes = 3;
        let phmm1 = dbg.to_phmm(param1);

        // 2
        let mut param2 = dataset.params();
        param2.n_active_nodes = 200;
        let phmm2 = dbg.to_phmm(param2);

        // Test score difference of these 3 phmms for each read..
        //
        for read in dataset.reads() {
            // o0: use max-ratio
            let (o0, t0) = timer(|| phmm0.run_sparse_adaptive(read, true));
            let p0f = o0.to_full_prob_forward();
            let p0b = o0.to_full_prob_backward();
            // o1: bad example, n_active_nodes is too small
            let (o1, t1) = timer(|| phmm1.run_sparse_adaptive(read, false));
            let p1f = o1.to_full_prob_forward();
            let p1b = o1.to_full_prob_backward();
            // o2: true score, n_active_nodes is large
            let (o2, t2) = timer(|| phmm2.run_sparse_adaptive(read, false));
            let p2f = o2.to_full_prob_forward();
            let p2b = o2.to_full_prob_backward();

            // Observations
            // * p1 can be wrong (insufficient n_active_nodes)
            // * p0 should be equal to p2 (max-ratio is accurate!)
            // * p0 should be faster than p2 (because of adaptive depth)
            println!("0 o={} {} t={}", p0f, p0b, t0);
            println!("1 o={} {} t={}", p1f, p1b, t1);
            println!("2 o={} {} t={}", p2f, p2b, t2);

            // Check 1:
            // p0 (use max ratio) vs p2 (ground truth)
            let diff = p0f.log_diff(p2f);
            assert!(diff < 0.0001);

            // Check 2:
            // p0 forward vs backward
            let diff = p0f.log_diff(p0b);
            assert!(diff < 0.0001);
        }
    }

    ///
    /// compare n_warmup adaptive or not
    ///
    fn test_read_dbg_n_warmup(genome: Genome, k: usize) {
        // create read DBG
        let (dataset, dbg) = dataset(genome, k);

        // phmm0: current fixed n_warmup
        let mut param0 = dataset.params();
        param0.n_warmup = k;
        param0.n_active_nodes = 200;
        let phmm0 = dbg.to_phmm(param0);

        // phmm1: adaptive n_warmup
        let mut param1 = dataset.params();
        param1.n_warmup = k;
        let phmm1 = dbg.to_phmm(param1);

        for read in dataset.reads() {
            let (o0, t0) = timer(|| phmm0.run_sparse_adaptive(read, false));
            let (o1, t1) = timer(|| phmm1.run_sparse_adaptive(read, true));
            println!(
                "0 pf={} pb={} t={}",
                o0.to_full_prob_forward(),
                o0.to_full_prob_backward(),
                t0,
            );
            println!(
                "1 pf={} pb={} t={}",
                o1.to_full_prob_forward(),
                o1.to_full_prob_backward(),
                t1,
            );

            // Check 1: o0 vs o1 forward
            let diff_01 = o0
                .to_full_prob_forward()
                .log_diff(o1.to_full_prob_forward());
            assert!(diff_01 < 0.0001);

            // Check 2: forward vs backward
            let diff_fb = o1
                .to_full_prob_forward()
                .log_diff(o1.to_full_prob_backward());
            assert!(diff_fb < 0.01);
            println!("diff 01={} fb={}", diff_01, diff_fb);
        }
    }

    #[test]
    fn u1k_read_dbg_mapping() {
        test_read_dbg_mapping(u1k(), 40)
    }
    #[test]
    fn u100_read_dbg_mapping() {
        test_read_dbg_mapping(u100(), 40)
    }
    #[test]
    fn u20_read_dbg_mapping() {
        test_read_dbg_mapping(u20(), 40)
    }
    #[test]
    fn u20_read_dbg_n_warmup() {
        println!("k=40");
        test_read_dbg_n_warmup(u20(), 40);
        println!("k=100");
        test_read_dbg_n_warmup(u20(), 100);
    }
    #[test]
    fn u1k_read_dbg_n_warmup() {
        println!("k=40");
        test_read_dbg_n_warmup(u1k(), 40);
        println!("k=100");
        test_read_dbg_n_warmup(u1k(), 100);

        // TODO
        // test k=500 or k=1000 but read dbg is inaccurate.
    }
    #[test]
    fn u20n200_read_dbg_max_ratio() {
        test_read_dbg_max_ratio(u20n200(), 40);
    }
    #[test]
    fn u20n200_read_dbg_n_warmup() {
        test_read_dbg_n_warmup(u20n200(), 40);
    }
}
