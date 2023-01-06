use dbgphmm::inspect::generate_500bp_case;
use std::fs::File;
use std::io::prelude::*;
fn main() {
    let (dataset, dbg_true, dbg_opt) = generate_500bp_case();
    dataset.show_reads_with_genome();
    println!("{}", dbg_true);
    println!("{}", dbg_opt);
    let copy_nums_true = dbg_true.to_node_copy_nums();
    let copy_nums_opt = dbg_opt.to_node_copy_nums();

    let json = dbg_true.to_cytoscape_with_info(
        |node| {
            Some(format!(
                "{}x {}x",
                copy_nums_true[node], copy_nums_opt[node]
            ))
        },
        None,
    );
    let mut file = File::create("inspect500bp_true_k12.json").unwrap();
    writeln!(file, "{}", json).unwrap();
}
