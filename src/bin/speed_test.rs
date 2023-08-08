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

    // create read DBG
    // dbg.to_gfa_file("read_dbg.gfa");
    let (mapping, t) = timer(|| dbg.generate_mappings(dataset.params(), dataset.reads(), None));
    println!("mapping created in t={}", t);
    for (r, m) in mapping.into_iter().enumerate() {
        for i in 0..m.len() {
            let l = m.nodes(i).len();
            println!("r={} i={} l={}", r, i, l);
        }
    }

    return;

    let phmm = dbg.to_phmm(dataset.params());

    //
    // (1) full-prob is same with/with-out mapping
    // likelihood of reads
    //
    // 0: without mapping, score only
    let (p0, t0) = timer(|| dbg.to_likelihood(dataset.params(), dataset.reads(), None));
    // 1: with mapping, score only
    let (p1, t1) = timer(|| dbg.to_likelihood(dataset.params(), dataset.reads(), Some(&mapping)));
    // 2: without mapping and do not use score only calculation
    let (p2, t2) = timer(|| phmm.to_full_prob_sparse(dataset.reads(), false));
    // 3:
    let (p3, t3) = timer(|| phmm.to_full_prob_sparse(dataset.reads(), true));
    let diff0 = p0.log_diff(p2);
    let diff1 = p1.log_diff(p2);
    println!(
        "p0={} p1={} p2={} p3={} diff0={} diff1={} t0={} t1={} t2={} t3={}",
        p0, p1, p2, p3, diff0, diff1, t0, t1, t2, t3
    );
}
