//!
//! Profile HMM speed test using large tandem repeat DBG
//!
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
    use crate::e2e::{generate_dataset, ReadType};
    use crate::genome;
    use crate::multi_dbg::MultiDbg;
    use crate::prelude::*;
    use crate::prob::{lp, p};
    use crate::utils::timer;
    use std::io::prelude::*;

    /// unit 1kb x 1000 times = 1MB diploid
    fn g1m() -> (Genome, usize) {
        let (genome, genome_size) = genome::tandem_repeat_polyploid_with_unique_homo_ends(
            1_000, 1_000, 0, 1_000, 2, 0.01, 0,
        );
        (genome, genome_size)
    }
    /// unit 20b x 50 times = 1kB diploid
    fn g1k() -> (Genome, usize) {
        let (genome, genome_size) =
            genome::tandem_repeat_polyploid_with_unique_homo_ends(20, 50, 0, 50, 2, 0.01, 0);
        (genome, genome_size)
    }

    #[test]
    #[ignore]
    fn test_g1m() {
        let (genome, genome_size) = g1m();
        let k = 40;
        let (dbg, time): (SimpleDbg<VecKmer>, _) =
            timer(|| SimpleDbg::from_styled_seqs(k, &genome));
        // ~0.5s in M1 mac
        println!(
            "dbg created in {}ms: n={} e={}",
            time,
            dbg.n_nodes(),
            dbg.n_edges()
        );
        assert_eq!(dbg.n_nodes(), 204734);
        let param = PHMMParams::uniform(0.01);
        let (phmm, time) = timer(|| dbg.to_phmm(param));
        // ~0.01s in M1 mac
        println!("phmm converted in {}ms", time);

        let (cedbg, time) = timer(|| dbg.to_compact_edbg_graph());
        println!("compact edbg converted in {}ms", time);

        // ~2s
        let (mdbg, time): (MultiDbg, _) = timer(|| dbg.into());
        println!("mdbg converted in {}ms", time);
        // mdbg.to_dbg_file("g1m.dbg");

        //
        // test using full-length error-free read
        //
        let (p, time) = timer(|| phmm.to_full_prob_sparse(&genome));
        // ~368sec (10min) in release on m1_mac
        println!("p={} t={}", p, time);
        let p_true = lp(-105736.2);
        assert!(p.log_diff(p_true) < 1.0);
    }

    #[test]
    fn test_g1k() {
        let (genome, genome_size) = g1k();
        let k = 40;
        let (dbg, time): (SimpleDbg<VecKmer>, _) =
            timer(|| SimpleDbg::from_styled_seqs(k, &genome));
        // ~9ms in M1 mac
        println!(
            "dbg created in {}ms: n={} e={}",
            time,
            dbg.n_nodes(),
            dbg.n_edges()
        );
        assert_eq!(dbg.n_nodes(), 512);
        let param = PHMMParams::uniform(0.01);
        let (phmm, time) = timer(|| dbg.to_phmm(param));
        // ~0ms in M1 mac
        println!("phmm converted in {}ms", time);

        // fragment reads
        let dataset = generate_dataset(
            genome.clone(),
            genome_size,
            0,
            2, // 2x
            500,
            ReadType::FragmentWithRevComp,
            param,
        );
        dataset.show_reads_with_genome();

        //
        // test using full-length error-free read
        //
        let p_true = lp(-138.0);
        // dense ~0.6sec
        let (p, time) = timer(|| phmm.to_full_prob(&genome));
        println!("p={} t={}", p, time);
        assert!(p.log_diff(p_true) < 1.0);
        // sparse ~0.3sec
        let (p, time) = timer(|| phmm.to_full_prob_sparse(&genome));
        println!("p={} t={}", p, time);
        assert!(p.log_diff(p_true) < 1.0);
        let (p, time) = timer(|| phmm.to_full_prob_sparse_backward(&genome));
        println!("p={} t={}", p, time);
        assert!(p.log_diff(p_true) < 1.0);

        //
        // test using fragment and errorneous read
        //
        println!("read");
        // dense
        let (p, time) = timer(|| phmm.to_full_prob(dataset.reads()));
        let p_true = lp(-1448.0);
        println!("p={} t={}", p, time);
        assert!(p.log_diff(p_true) < 1.0);
        // sparse
        let (p, time) = timer(|| phmm.to_full_prob_sparse(dataset.reads()));
        let p_true = lp(-1448.0);
        println!("p={} t={}", p, time);
        assert!(p.log_diff(p_true) < 1.0);

        let mut file = std::fs::File::create("g1k_dense_vs_sparse_for_reads.tsv").unwrap();
        writeln!(
            file,
            "# read index base origin log_diff_end p_end_dense p_end_sparse"
        )
        .unwrap();

        for (r, read) in dataset.reads().iter().enumerate() {
            // dense
            let dense = phmm.forward(read);
            let sparse = phmm.forward_sparse(read);
            let n = read.len();
            for i in 0..n {
                // println!("dense");
                // println!("{}", dense.table(i));
                // println!("sparse");
                // println!("{}", sparse.table(i));
                let td = dense.table(i);
                let ts = sparse.table(i);
                writeln!(
                    file,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    r,
                    i,
                    read.seq()[i] as char,
                    read.origins()[i],
                    td.e.log_diff(ts.e),
                    td.e,
                    ts.e,
                    td.to_summary_string_n(5, |_| String::new()),
                    ts.to_summary_string_n(5, |_| String::new()),
                    td.diff(ts),
                )
                .unwrap();
            }
        }

        // backward
        // let (p, time) = timer(|| phmm.to_full_prob_sparse_backward(dataset.reads()));
        // println!("p={} t={}", p, time);
        // assert!(p.log_diff(p_true) < 1.0);
    }
}
