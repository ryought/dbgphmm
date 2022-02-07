use dbgphmm::graph::genome_graph;

fn main() {
    let g = genome_graph::mock();
    // println!("{}", g);
    let sg = g.to_seq_graph();
    println!("{}", sg);
}
