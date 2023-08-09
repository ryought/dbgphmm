use dbgphmm::{
    e2e::{generate_dataset, ReadType},
    genome,
    hmmv2::params::PHMMParams,
    multi_dbg::MultiDbg,
    utils::{check_memory_usage, timer},
};

fn main() {
    let genome =
        genome::tandem_repeat_polyploid_with_unique_homo_ends(20, 200, 0, 0.02, 0, 300, 2, 0.02, 0);
    let dataset = generate_dataset(
        genome,
        0,
        20,                            // 20x coverage
        1000,                          // 1000bp reads
        ReadType::FragmentWithRevComp, // fragment read
        PHMMParams::uniform(0.001),    // 0.1% HiFi error
    );
    let (dbg, t) = timer(|| MultiDbg::create_draft_from_dataset(40, &dataset));

    let phmm = dbg.to_phmm(dataset.params());

    /*
    for (i, read) in dataset.reads().into_iter().enumerate() {
        let (p, t) = timer(|| phmm.forward_sparse(read, true));
        for j in 0..p.n_emissions() {
            println!(
                "j={j} n_top_nodes={} is_dense={}",
                p.table(j).top_nodes_by_score_ratio(30.0).len(),
                p.table(j).is_dense()
            );
        }
    }
    */

    // create read DBG
    // dbg.to_gfa_file("read_dbg.gfa");
    let (mappings, t) = timer(|| dbg.generate_mappings(dataset.params(), dataset.reads(), None));
    println!("mapping created in t={}", t);
    for (r, m) in mappings.into_iter().enumerate() {
        for i in 0..m.len() {
            let l = m.nodes(i).len();
            println!("r={} i={} l={}", r, i, l);
        }
    }

    /*
    //
    // (1) full-prob is same with/with-out mapping
    // likelihood of reads
    //
    // 0: without mapping, score only
    let (p0, t0) = timer(|| dbg.to_likelihood(dataset.params(), dataset.reads(), None));
    // 1: with mapping, score only
    let (p1, t1) = timer(|| dbg.to_likelihood(dataset.params(), dataset.reads(), Some(&mappings)));
    // 2: without mapping and do not use score only calculation
    let (p2, t2) = timer(|| phmm.to_full_prob_sparse(dataset.reads(), false));
    // 3:
    let (p3, t3) = timer(|| phmm.to_full_prob_sparse(dataset.reads(), true));
    let diff0 = p0.log_diff(p2);
    let diff1 = p1.log_diff(p2);
    println!(
        "p0={} p1={} p2={} p3={} diff0={} diff1={}",
        p0, p1, p2, p3, diff0, diff1
    );
    println!("t0={t0} t1={t1} t2={t2} t3={t3}");
    */

    //
    // (2)
    //
    for (i, read) in dataset.reads().into_iter().enumerate() {
        let (p1, t1) = timer(|| phmm.forward_with_mapping_score_only(read, &mappings[i]));
        let (p2, t2) = timer(|| phmm.forward_sparse_score_only(read, true));
        let (p3, t3) = timer(|| phmm.forward_sparse_score_only(read, false));
        let (p4, t4) = timer(|| phmm.forward_sparse_v0(read, true));
        let (p5, t5) = timer(|| phmm.forward_sparse_v0(read, false));
        let (p6, t6) = timer(|| phmm.forward_sparse(read, true));
        let (p7, t7) = timer(|| phmm.forward_sparse(read, false));

        println!("{t1}\t{t2}\t{t3}\t{t4}\t{t5}\t{t6}\t{t7}");
    }
}
