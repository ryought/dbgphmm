//!
//! Markov chain monte carlo methods
//!
//! Random sample
//!
use crate::dbg::dbg::{Dbg, DbgEdge, DbgNode, EdgeCopyNums, NodeCopyNums};
use crate::dbg::edge_centric::{EDbg, EDbgEdge, EDbgNode};
use crate::graph::shortest_cycle::shortest_cycle;
use crate::min_flow::residue::{
    change_flow_along_edges, ResidueDirection, ResidueEdge, ResidueGraph,
};
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
use petgraph::prelude::*;
use petgraph_algos::common::nodes_to_edges;

///
/// get a disturbed neighboring copy numbers
///
fn get_neighbor_copy_nums<N: DbgNode, E: DbgEdge>(
    dbg: &Dbg<N, E>,
) -> Vec<(NodeIndex, ResidueDirection, NodeCopyNums)> {
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
                    NodeIndex::new(target.index()), // assuming edge in edbg == node in dbg
                    direction,
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
fn get_reverse_edge(graph: &ResidueGraph<usize>, edge: EdgeIndex) -> Option<EdgeIndex> {
    let weight = graph.edge_weight(edge).unwrap();
    let target = weight.target;
    let direction = weight.direction;

    graph.edge_indices().find(|&e| {
        let ew = graph.edge_weight(e).unwrap();
        ew.target == target && ew.direction != direction
    })
}

fn edbg_to_residue<N: EDbgNode, E: EDbgEdge>(edbg: &EDbg<N, E>) -> ResidueGraph<usize> {
    let mut rg = ResidueGraph::new();
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
    rg: &ResidueGraph<usize>,
    cycle: &[EdgeIndex],
) -> NodeCopyNums {
    let flow: EdgeCopyNums = copy_nums.clone().switch_index();
    let new_flow = change_flow_along_edges(&flow.into(), rg, cycle, 1);
    let copy_nums: EdgeCopyNums = new_flow.into();
    copy_nums.switch_index()
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
        let copy_nums_list = get_neighbor_copy_nums(&dbg);

        // (0) the number of candidate copy_nums is correct?
        assert_eq!(copy_nums_list.len(), dbg.n_nodes() * 2);

        for (v, dir, copy_nums) in copy_nums_list {
            println!("v={:?} direction={} copy_nums={}", v, dir, copy_nums);

            // (1) copy_nums have actually increased/decreased the copy_num of the targeted node?
            match dir {
                ResidueDirection::Up => {
                    assert_eq!(copy_nums[v], 2);
                }
                ResidueDirection::Down => {
                    assert_eq!(copy_nums[v], 0);
                }
            };

            // TODO to debug
            dbg.draw_with_vecs(&[&copy_nums], &[]);

            // (2) the resulting copy_nums is valid?
            let mut new_dbg = dbg.clone();
            new_dbg.set_node_copy_nums(&copy_nums);
            let is_valid = new_dbg.has_consistent_node_copy_nums();
            assert!(is_valid);
        }
    }
}
