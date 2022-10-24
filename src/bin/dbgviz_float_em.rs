use dbgphmm::common::{ni, Reads};
use dbgphmm::dbg::float::{q_score_diff_exact, CopyDensity, FloatDbg, FloatDbgEdge, FloatDbgNode};
use dbgphmm::dbg::mocks::mock_intersection_small;
use dbgphmm::e2e::{generate_dataset, Dataset, ReadType};
use dbgphmm::em::float::{
    em, em_with_upgrade, inspect_density_histgram, inspect_freqs_histgram, run, shrink_nodes,
};
use dbgphmm::genome;
use dbgphmm::io::write_string;
use dbgphmm::prelude::*;
use dbgphmm::vector::{DenseStorage, NodeVec};
use rayon::prelude::*;

fn run_small() {
    let dbg = mock_intersection_small();
    let genome_size = dbg.genome_size() as CopyDensity;
    let reads = Reads::from(vec![b"ATAGCT".to_vec()]);
    let fdbg = FloatDbg::from_dbg(&dbg);
    let params = PHMMParams::zero_error();

    let result = em(&fdbg, &reads, &params, genome_size, 0.1, 10, 10);
    let historys = result.to_node_historys();

    let json = dbg.to_cytoscape_with_attrs_and_historys(&[], &[], &historys);
    // let (edge_freqs, init_freqs, p) = e_step(&fdbg, &reads, &params);
    println!("{}", json);
}

fn run_simple() {
    // let (genome, genome_size) = genome::simple(100, 5);
    let (genome, genome_size) = genome::tandem_repeat_haploid(50, 4, 0.05, 0, 0);
    eprintln!("{}", genome[0]);
    let coverage = 10;
    let param = PHMMParams::uniform(0.001);
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

    // dbg_raw
    let dbg_raw = dataset.dbg_raw;

    // dbg_true
    let mut dbg_true = dbg_raw.clone();
    let (copy_nums_true, _) = dbg_true.to_copy_nums_of_styled_seqs(&genome).unwrap();
    dbg_true.set_node_copy_nums(&copy_nums_true);

    let run_em = true;
    let run_shrink = true;

    if run_em {
        eprintln!("{}", dbg_raw.to_node_copy_nums());
        let mut fdbg = FloatDbg::from_dbg(&dbg_raw);
        eprintln!("n_nodes={}", fdbg.n_nodes());
        eprintln!("n_edges={}", fdbg.n_edges());
        fdbg.scale_density(genome_size as CopyDensity / dbg_raw.genome_size() as CopyDensity);
        eprintln!("{}", fdbg.to_node_copy_densities());

        let result = em(
            &fdbg,
            &dataset.reads,
            &param,
            genome_size as CopyDensity,
            0.001,
            50,
            50,
        );
        let mut historys = result.to_node_historys();
        let final_fdbg = &result.to_final_dbg();
        result.inspect_stop_reason();

        if let Some(final_fdbg) = final_fdbg {
            // inspect histogram
            final_fdbg.inspect_freqs_histogram(&genome);

            let shrink_min_density = 0.05;
            eprintln!("n_red={}", fdbg.n_redundant_nodes(shrink_min_density));
            eprintln!("n_dead={}", fdbg.n_dead_nodes());

            // shrink
            if run_shrink {
                let shrinked_densities = shrink_nodes(&final_fdbg, shrink_min_density);
                let mut fdbg_shrinked = fdbg.clone();
                fdbg_shrinked.set_node_copy_densities(&shrinked_densities);
                eprintln!(
                    "n_red={}",
                    fdbg_shrinked.n_redundant_nodes(shrink_min_density)
                );
                eprintln!("n_dead={}", fdbg_shrinked.n_dead_nodes());

                historys.push((format!("shrinked"), shrinked_densities))
            }
        }

        let json = dbg_true.to_cytoscape_with_attrs_and_historys(&[], &[], &historys);
        println!("{}", json);
    } else {
        let phmm = dbg_raw.to_phmm(param.clone());
        let (nf, p) = phmm.to_node_freqs_parallel(&dataset.reads);
        inspect_freqs_histgram(&dbg_true, &nf);
        let json = dbg_true.to_cytoscape_with_attrs_and_historys(
            &[],
            &[],
            &[("node_freqs".to_string(), nf)],
        );
        println!("{}", json);
    }
}

fn run_mcmc() {
    // let (genome, genome_size) = genome::simple(100, 5);
    let (genome, genome_size) = genome::tandem_repeat_haploid(50, 4, 0.05, 0, 0);
    let coverage = 10;
    let param = PHMMParams::uniform(0.001);
    let dataset = generate_dataset(
        genome.clone(),
        genome_size,
        0,
        param.clone(),
        coverage,
        2000,
        ReadType::FullLength,
        8,
        40,
    );
    let dbg_raw = dataset.dbg_raw.clone();
    let mut dbg_true = dbg_raw.clone();
    let (copy_nums_true, _) = dbg_true.to_copy_nums_of_styled_seqs(&genome).unwrap();
    dbg_true.set_node_copy_nums(&copy_nums_true);

    let mut neighbors = dbg_true.neighbor_copy_nums();
    neighbors.push(copy_nums_true);
    neighbors.par_iter().for_each(|copy_nums| {
        let mut dbg = dbg_true.clone();
        dbg.set_node_copy_nums(copy_nums);
        let p = dbg.to_full_prob(param, &dataset.reads);
        println!(
            "{}\t{}\t{}\t{}",
            dbg.genome_size(),
            p.to_log_value(),
            copy_nums,
            dbg,
        );
    })
}

fn run_upgrade() {
    // let (genome, genome_size) = genome::simple(100, 5);
    let (genome, genome_size) = genome::tandem_repeat_haploid(50, 4, 0.05, 0, 0);
    eprintln!("{}", genome[0]);
    let coverage = 20;
    let param = PHMMParams::uniform(0.01);
    let dataset = generate_dataset(
        genome.clone(),
        genome_size,
        0,
        param.clone(),
        coverage,
        2000,
        ReadType::FullLength,
        8,
        64,
    );
    // dataset.show_reads();

    // dbg_raw
    let dbg_raw = dataset.dbg_raw;

    // dbg_true
    let mut dbg_true = dbg_raw.clone();
    let (copy_nums_true, _) = dbg_true.to_copy_nums_of_styled_seqs(&genome).unwrap();
    dbg_true.set_node_copy_nums(&copy_nums_true);

    let mut fdbg = FloatDbg::from_dbg(&dbg_raw);
    eprintln!("n_nodes={}", fdbg.n_nodes());
    eprintln!("n_edges={}", fdbg.n_edges());
    fdbg.scale_by_total_density(genome_size as CopyDensity);

    run(
        &fdbg,
        &dataset.reads,
        &param,
        genome_size as CopyDensity,
        0.001,
        50,
        50,
        0.1,
        40,
        |((init, p_init), (opt, p_opt), (shrinked, p_shrinked))| {
            init.benchmark(&genome, *p_init);
            opt.benchmark(&genome, *p_opt);
            shrinked.benchmark(&genome, *p_shrinked);
        },
    );

    // let json = dbg_true.to_cytoscape_with_attrs_and_historys(&[], &[], &historys);
    // println!("{}", json);
}

fn run_frag() {
    // let (genome, genome_size) = genome::simple_circular(100, 5);
    let (genome, genome_size) = genome::simple(100, 5);
    // let (genome, genome_size) = genome::tandem_repeat_haploid(50, 4, 0.05, 0, 0);
    eprintln!("{}", genome[0]);
    let coverage = 20;
    let param = PHMMParams::uniform(0.001);
    let dataset = generate_dataset(
        genome.clone(),
        genome_size,
        0,
        param.clone(),
        coverage,
        50,
        ReadType::FixedSizeFragment,
        8,
        64,
    );
    dataset.show_reads();

    let json = dataset.dbg_raw.to_cytoscape();
    write_string("linear_raw.json", &json);

    // optimize
    let mut fdbg = FloatDbg::from_dbg(&dataset.dbg_raw);
    eprintln!("n_nodes={}", fdbg.n_nodes());
    eprintln!("n_edges={}", fdbg.n_edges());
    fdbg.scale_by_total_density(genome_size as CopyDensity);

    let upgrade = false;

    // let phmm = dataset.dbg_raw.to_phmm(param);
    // let output = phmm.run(&dataset.reads[0]);
    // dataset
    //     .dbg_raw
    //     .show_mapping_summary(&dataset.reads[0], &output);

    if upgrade {
        // upgrade
        run(
            &fdbg,
            &dataset.reads,
            &param,
            genome_size as CopyDensity,
            0.001,
            50,
            50,
            0.05,
            40,
            |((init, p_init), (opt, p_opt), (shrinked, p_shrinked))| {
                init.benchmark(&genome, *p_init);
                opt.benchmark(&genome, *p_opt);
                shrinked.benchmark(&genome, *p_shrinked);
            },
        );
    } else {
        // single
        let result = em(
            &fdbg,
            &dataset.reads,
            &param,
            genome_size as CopyDensity,
            0.001,
            50,
            50,
        );
        let mut historys = result.to_node_historys();

        let mut dbg_true = dataset.dbg_raw.clone();
        let (copy_nums_true, _) = dbg_true.to_copy_nums_of_styled_seqs(&genome).unwrap();
        dbg_true.set_node_copy_nums(&copy_nums_true);
        let json = dbg_true.to_cytoscape_with_attrs_and_historys(&[], &[], &historys);
        write_string("linear_em.json", &json);
    }
}

fn main() {
    run_mcmc();
    // run_frag();
    // run_upgrade();
    // run_simple();
}
