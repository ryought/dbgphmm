use dbgphmm::graph::mocks::mock;

fn main() {
    let g = mock();
    // println!("{}", g);
    let sg = g.to_seq_graph();
    let phmm = sg.to_phmm();
    println!("{}", phmm);
}
