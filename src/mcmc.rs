//!
//! Markov chain monte carlo methods
//!
//! Random sample
//!
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, NodeCopyNums};
use crate::dbg::edge_centric::{EDbg, EDbgEdge, EDbgNode};
use crate::graph::float_weight::nodes_to_edges;
use crate::graph::shortest_cycle::shortest_cycle;
use crate::min_flow::residue::{
    change_flow_along_edges, ResidueDirection, ResidueEdge, ResidueGraph,
};
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
use petgraph::prelude::*;

///
/// get a disturbed neighboring copy numbers
///
fn get_neighbor_copy_nums<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
) -> Vec<(NodeIndex, bool, NodeCopyNums)> {
    // a list of copy_nums
    let mut ret = Vec::new();

    let edbg = dbg.to_edbg();
    let copy_nums = dbg.to_node_copy_nums();
    let rg = edbg_to_residue(&edbg);

    for e_for in rg.edge_indices() {
        let ew = rg.edge_weight(e_for).unwrap();
        let target = ew.target;
        let direction = ew.direction;

        let e_rev = get_reverse_edge(&rg, e_for);
        let cycle = shortest_cycle(&rg, e_for, e_rev);
        match cycle {
            Some(cycle) => {
                let cycle = nodes_to_edges(&rg, &cycle, |graph, v, w| {
                    graph.edges_connecting(v, w).next().unwrap().id()
                });
                let new_copy_nums = apply_residual_cycles_to_copy_nums(&copy_nums, &rg, &cycle);

                ret.push((
                    NodeIndex::new(target.index()), // XXX edge in edbg == node in dbg
                    direction == ResidueDirection::Down,
                    new_copy_nums,
                ))
            }
            None => {}
        }
    }

    ret
}

///
/// get a reverse edge, which has the same `target` and the opposite `direction`.
///
fn get_reverse_edge(graph: &ResidueGraph, edge: EdgeIndex) -> Option<EdgeIndex> {
    let weight = graph.edge_weight(edge).unwrap();
    let target = weight.target;
    let direction = weight.direction;

    graph.edge_indices().find(|&e| {
        let ew = graph.edge_weight(e).unwrap();
        ew.target == target && ew.direction != direction
    })
}

fn edbg_to_residue<N: EDbgNode, E: EDbgEdge>(edbg: &EDbg<N, E>) -> ResidueGraph {
    let mut rg: ResidueGraph = ResidueGraph::new();
    // create two edges (Up and Down) for each edge
    for (e, v, w, ew) in edbg.edges() {
        let mut edges = Vec::new();
        let copy_num = ew.copy_num();
        if copy_num < usize::MAX {
            // up movable
            edges.push((v, w, ResidueEdge::new(1, 0.0, e, ResidueDirection::Up)));
        }
        if copy_num > 0 {
            // down movable
            edges.push((w, v, ResidueEdge::new(1, 0.0, e, ResidueDirection::Down)));
        }
        rg.extend_with_edges(&edges);
    }
    rg
}

///
///
///
fn apply_residual_cycles_to_copy_nums(
    copy_nums: &NodeCopyNums,
    rg: &ResidueGraph,
    cycle: &[EdgeIndex],
) -> NodeCopyNums {
    let flow = copy_nums.clone().switch_index();
    let new_flow = change_flow_along_edges(&flow, rg, cycle, 1);
    new_flow.switch_index()
}

trait Model {
    fn score(&self) -> f64;
}

struct Runner {
    // proposal: P,
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dbg::mocks::mock_intersection_small;

    #[test]
    fn dbg_neighbor_copy_nums() {
        let dbg = mock_intersection_small();
        println!("{}", dbg.to_dot());
        let copy_nums_list = get_neighbor_copy_nums(&dbg);
        for (v, is_down, copy_nums) in copy_nums_list {
            println!("v={:?} is_down={} copy_nums={}", v, is_down, copy_nums);
        }
    }
}
