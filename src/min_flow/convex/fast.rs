//!
//! Fast min-flow solver for convex cost functions
//!
//! # (1) naive
//!
//! * convert to the constant cost graph, by duplicating edges which has a constant cost of `f(i+1)-f(i)`.
//! * solve the graph by normal min-flow solver.
//!
//! # (2) fast
//!
//! * to find the initial valid flow, cost will be unit and...
//! * to find the min flow from the init flow, residual...
//!
use super::super::flow::{Flow, FlowEdge};
use super::ConvexCost;
use petgraph::graph::DiGraph;

/// create a new improved flow from current flow
/// by upgrading along the negative weight cycle in the residual graph
pub fn improve_flow_convex<N, E>(graph: &DiGraph<N, E>, flow: &Flow) -> Option<Flow>
where
    E: FlowEdge + ConvexCost,
{
    unimplemented!();
    /*
    let rg = flow_to_residue(graph, flow);

    // find negative weight cycles
    let path = find_negative_cycle_in_whole_graph(&rg);
    draw(&rg);

    match path {
        Some(nodes) => {
            let edges = node_list_to_edge_list(&rg, &nodes);

            // check if this is actually negative cycle
            assert!(is_negative_cycle(&rg, &edges));

            // apply these changes along the cycle to current flow
            let new_flow = apply_residual_edges_to_flow(&flow, &rg, &edges);
            println!("{:?}", new_flow);
            Some(new_flow)
        }
        None => None,
    }
    */
}
