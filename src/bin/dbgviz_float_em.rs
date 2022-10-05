use dbgphmm::common::{ni, Reads};
use dbgphmm::dbg::float::{q_score_diff_exact, CopyDensity, FloatDbg, FloatDbgEdge, FloatDbgNode};
use dbgphmm::dbg::mocks::mock_intersection_small;
use dbgphmm::e2e::{generate_dataset, Dataset, ReadType};
use dbgphmm::em::float::{
    em, em_result_to_final_dbg, em_result_to_node_historys, inspect_density_histgram,
};
use dbgphmm::genome;
use dbgphmm::prelude::*;
use dbgphmm::vector::{DenseStorage, NodeVec};

fn run_small() {
    let dbg = mock_intersection_small();
    let genome_size = dbg.genome_size() as CopyDensity;
    let reads = Reads::from(vec![b"ATAGCT".to_vec()]);
    let fdbg = FloatDbg::from_dbg(&dbg);
    let params = PHMMParams::zero_error();

    let result = em(&fdbg, &reads, &params, genome_size, 0.1, 10, 10);
    let historys = em_result_to_node_historys(&result);

    let json = dbg.to_cytoscape_with_attrs_and_historys(&[], &[], &historys);
    // let (edge_freqs, init_freqs, p) = e_step(&fdbg, &reads, &params);
    println!("{}", json);
}

fn run_simple() {
    // let (genome, genome_size) = genome::simple(100, 5);
    let (genome, genome_size) = genome::tandem_repeat_haploid(50, 4, 0.05, 0, 0);
    eprintln!("{}", genome[0]);
    let coverage = 10;
    let param = PHMMParams::uniform(0.01);
    let dataset = generate_dataset(
        genome.clone(),
        genome_size,
        0,
        param.clone(),
        coverage,
        2000,
        ReadType::FullLength,
        12,
        64,
    );
    // dataset.show_reads();

    let dbg_raw = dataset.dbg_raw;
    // let phmm = dbg_raw.to_phmm(param.clone());
    // let (nf, p) = phmm.to_node_freqs_parallel(&dataset.reads);
    let mut fdbg = FloatDbg::from_dbg(&dbg_raw);
    fdbg.scale_density(genome_size as CopyDensity / dbg_raw.genome_size() as CopyDensity);

    let result = em(
        &fdbg,
        &dataset.reads,
        &param,
        genome_size as CopyDensity,
        0.05,
        10,
        10,
    );
    let historys = em_result_to_node_historys(&result);
    let final_fdbg = em_result_to_final_dbg(&result);

    let mut dbg_true = dbg_raw.clone();
    let (copy_nums_true, _) = dbg_true.to_copy_nums_of_styled_seqs(&genome).unwrap();
    dbg_true.set_node_copy_nums(&copy_nums_true);

    if let Some(final_fdbg) = final_fdbg {
        inspect_density_histgram(&dbg_true, &final_fdbg);
    }

    // let json =
    //     dbg_true.to_cytoscape_with_attrs_and_historys(&[], &[], &[("node_freqs".to_string(), nf)]);
    let json = dbg_true.to_cytoscape_with_attrs_and_historys(&[], &[], &historys);
    println!("{}", json);
}

fn main() {
    run_simple();
}
