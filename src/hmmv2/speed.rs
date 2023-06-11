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

    /// unit 1kb x 1000 times = 1MB diploid
    fn g1m() -> Genome {
        genome::tandem_repeat_polyploid_with_unique_homo_ends(
            1_000, 1_000, 0, 0.0, 0, 1_000, 2, 0.01, 0,
        )
    }
    /// unit 20b x 50 times = 1kB diploid
    fn g1k() -> Genome {
        genome::tandem_repeat_polyploid_with_unique_homo_ends(20, 50, 0, 0.0, 0, 50, 2, 0.01, 0)
    }
    /// unit 20b x 500 times = 1kB diploid
    fn g10k() -> Genome {
        genome::tandem_repeat_polyploid_with_unique_homo_ends(20, 500, 0, 0.0, 0, 50, 2, 0.01, 0)
    }
    /// unit 1kB x 1 times = 1kB diploid (unique sequence)
    fn gu1k() -> Genome {
        genome::tandem_repeat_polyploid_with_unique_homo_ends(1000, 1, 0, 0.0, 0, 50, 2, 0.01, 0)
    }

    /// Forward/backward table
    ///
    ///
    fn table_to_active_nodes_with_prob(
        table: &PHMMTable,
        p0: Prob,
    ) -> (Vec<(usize, f64)>, Vec<(usize, f64)>, Vec<(usize, f64)>) {
        // sum of probability exceeds p0

        // m
        let mut ms = Vec::new();
        let p_m_sum: Prob = table.m.iter().map(|(_, p)| p).sum();
        let mut n_m = 0;
        let mut p_m = Prob::zero();
        for (node, p) in table.m.iter().sorted_by_key(|(_, p)| *p).rev() {
            p_m += p / p_m_sum;
            if p_m > p0 {
                break;
            }
            ms.push((node.index(), (p / p_m_sum).to_value()));
            n_m += 1;
        }

        // i
        let mut is = Vec::new();
        let p_i_sum: Prob = table.i.iter().map(|(_, p)| p).sum();
        let mut n_i = 0;
        let mut p_i = Prob::zero();
        for (node, p) in table.i.iter().sorted_by_key(|(_, p)| *p).rev() {
            p_i += p / p_i_sum;
            if p_i > p0 {
                break;
            }
            is.push((node.index(), (p / p_i_sum).to_value()));
            n_i += 1;
        }

        // d
        let mut ds = Vec::new();
        let p_d_sum: Prob = table.d.iter().map(|(_, p)| p).sum();
        let mut n_d = 0;
        let mut p_d = Prob::zero();
        for (node, p) in table.d.iter().sorted_by_key(|(_, p)| *p).rev() {
            p_d += p / p_d_sum;
            if p_d > p0 {
                break;
            }
            ds.push((node.index(), (p / p_d_sum).to_value()));
            n_d += 1;
        }

        // (n_m, n_i, n_d)
        (ms, is, ds)
    }

    fn format_states(states: &[(State, Prob)]) -> String {
        states
            .iter()
            .map(|(state, prob)| format!("{}:{:.3}", state, prob.to_value()))
            .join(",")
    }

    fn inspect_forward_vs_state(dataset: &Dataset, phmm: &PModel, key: &str) {
        let p0 = Prob::from_prob(0.999);
        // forward
        let mut file = std::fs::File::create(format!("{}_forward_vs_state.tsv", key)).unwrap();
        for (r, read) in dataset.reads().iter().enumerate() {
            let (output, t) = timer(|| phmm.run(read));
            println!("{} {}", r, t);
            let n = read.len();
            for i in 0..n {
                writeln!(
                    file,
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    r,
                    i,
                    read.seq()[i] as char,
                    read.origins()[i],
                    format_states(&output.forward.table(i).to_states_normalize(Some(p0), true)),
                    format_states(
                        &output
                            .to_emit_probs(i + 1)
                            .to_states_normalize(Some(p0), false)
                    ),
                )
                .unwrap();
            }
        }
    }

    fn inspect_sparse_vs_dense(dataset: &Dataset, phmm: &PModel, key: &str) {
        let p0 = Prob::from_prob(0.999);
        // forward
        let mut file =
            std::fs::File::create(format!("{}_dense_vs_sparse_for_reads_forward.tsv", key))
                .unwrap();
        for (r, read) in dataset.reads().iter().enumerate() {
            // dense
            let (dense, t_dense) = timer(|| phmm.forward(read));
            let (sparse, t_sparse) = timer(|| phmm.forward_sparse(read, false));
            let (p, t_sparse_score) = timer(|| phmm.forward_sparse_score_only(read));
            println!("{} {} {}", t_dense, t_sparse, t_sparse_score);
            let n = read.len();
            for i in 0..n {
                // println!("dense");
                // println!("{}", dense.table(i));
                // println!("sparse");
                // println!("{}", sparse.table(i));
                let td = dense.table(i);
                let ts = sparse.table(i);
                let (n_m, n_i, n_d) = table_to_active_nodes_with_prob(&td, p0);
                writeln!(
                    file,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:?}\t{:?}\t{:?}\t{}\t{}",
                    r,
                    i,
                    read.seq()[i] as char,
                    read.origins()[i],
                    td.e.log_diff(ts.e),
                    td.e,
                    ts.e,
                    n_m,
                    n_i,
                    n_d,
                    td.to_summary_string_n(5, |_| String::new()),
                    ts.to_summary_string_n(5, |_| String::new()),
                )
                .unwrap();
            }
        }

        // backward
        let mut file =
            std::fs::File::create(format!("{}_dense_vs_sparse_for_reads_backward.tsv", key))
                .unwrap();
        for (r, read) in dataset.reads().iter().enumerate() {
            // dense
            let dense = phmm.backward(read);
            let sparse = phmm.backward_sparse(read);
            let n = read.len();
            for i in 0..n {
                let td = dense.table(i);
                let ts = sparse.table(i);
                let (n_m, n_i, n_d) = table_to_active_nodes_with_prob(&td, p0);
                writeln!(
                    file,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:?}\t{:?}\t{:?}\t{}\t{}",
                    r,
                    i,
                    read.seq()[i] as char,
                    read.origins()[i],
                    td.mb.log_diff(ts.mb),
                    td.mb,
                    ts.mb,
                    n_m,
                    n_i,
                    n_d,
                    td.to_summary_string_n(5, |_| String::new()),
                    ts.to_summary_string_n(5, |_| String::new()),
                )
                .unwrap();
            }
        }
    }

    #[test]
    #[ignore]
    fn test_g1m() {
        let genome = g1m();
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
        let mut param = PHMMParams::uniform(0.01);
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
        let p_true = lp(-105721.1);
        assert!(p.log_diff(p_true) < 1.0);
    }

    #[test]
    #[ignore]
    fn test_g10k() {
        let genome = g10k();
        let k = 40;
        let dbg: SimpleDbg<VecKmer> = SimpleDbg::from_styled_seqs(k, &genome);
        println!("{:?}", dbg.degree_stats());
        let param = PHMMParams::uniform(0.01);
        let phmm = dbg.to_phmm(param);

        let dataset = generate_dataset(
            genome,
            0,
            2, // 2x
            500,
            ReadType::FragmentWithRevComp,
            param,
        );
        inspect_sparse_vs_dense(&dataset, &phmm, "g10k");
    }

    #[test]
    #[ignore = "39sec"]
    fn test_g1k() {
        let genome = g1k();
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
        let (mut phmm, time) = timer(|| dbg.to_phmm(param));
        // ~0ms in M1 mac
        println!("phmm converted in {}ms", time);

        // fragment reads
        let dataset = generate_dataset(
            genome.clone(),
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
        // sparse ~0.3sec
        let (p, time) = timer(|| phmm.to_full_prob_sparse(&genome));
        let (p, time) = timer(|| phmm.to_full_prob_sparse_backward(&genome));
        println!("p={} t={}", p, time);
        assert!(p.log_diff(p_true) < 1.0);

        //
        // test using fragment and errorneous read
        //
        println!("read");
        // dense ~1.6s
        let (p, time) = timer(|| phmm.to_full_prob(dataset.reads()));
        let p_true = lp(-1448.0);
        println!("p={} t={}", p, time);
        assert!(p.log_diff(p_true) < 1.0);
        // sparse ~0.5s
        let (p, time) = timer(|| phmm.to_full_prob_sparse(dataset.reads()));
        let p_true = lp(-1448.0);
        println!("p={} t={}", p, time);
        assert!(p.log_diff(p_true) < 1.0);
        // sparse backward ~0.6s
        // n_active_nodes = 40 causes inaccurate prob
        phmm.param.n_active_nodes = 40;
        let (p, time) = timer(|| phmm.to_full_prob_sparse_backward(dataset.reads()));
        let p_true = lp(-1486.0);
        println!("p={} t={}", p, time);
        assert!(p.log_diff(p_true) < 1.0);
        // n = 100 is ok.
        phmm.param.n_active_nodes = 100;
        let (p, time) = timer(|| phmm.to_full_prob_sparse_backward(dataset.reads()));
        let p_true = lp(-1448.0);
        println!("p={} t={}", p, time);
        assert!(p.log_diff(p_true) < 1.0);

        // let (p, time) = timer(|| phmm.to_full_prob_sparse_backward(dataset.reads()));
        // println!("p={} t={}", p, time);
        // assert!(p.log_diff(p_true) < 1.0);

        inspect_sparse_vs_dense(&dataset, &phmm, "g1k");

        inspect_forward_vs_state(&dataset, &phmm, "g1k");
    }
}
