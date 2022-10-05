use dbgphmm::common::Reads;
use dbgphmm::dbg::mocks::mock_intersection;

fn main() {
    let mut dbg = mock_intersection();
    let json = dbg.to_cytoscape();
    println!("{}", json);
}
