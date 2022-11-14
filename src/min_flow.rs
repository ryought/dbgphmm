pub mod common;
pub mod convex;
pub mod flow;
pub mod mocks;
pub mod residue;
pub mod utils;
pub mod zero_demand;

use crate::graph::cycle::CycleWithDir;
use convex::{is_convex_cost_flow_graph, restore_convex_flow, to_fixed_flow_graph, ConvexCost};
pub use flow::total_cost;
use flow::{assert_valid_flow, is_valid_flow, ConstCost, Flow, FlowEdge, FlowGraphRaw};
use petgraph::graph::DiGraph;
use residue::{
    enumerate_neighboring_flows_in_residue, flow_to_residue_convex, improve_flow,
    improve_flow_convex, CycleDetectMethod, UpdateInfo,
};
use utils::draw_with_flow;
use zero_demand::{find_initial_flow, is_zero_demand_flow_graph};

/// type of a flow (on edges) in min-flow.
pub type FlowRate = usize;

use std::iter::{Step, Sum};
use std::ops::{Add, AddAssign, Div, Mul, Sub};
///
/// generic FlowRate
///
pub trait FlowRateLike:
    Copy
    + PartialEq
    + PartialOrd
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + AddAssign
    + Sum
    + Default
    + std::fmt::Debug
    + std::fmt::Display
{
    /// zero value = 0
    fn zero() -> Self;
    /// unit value = 1
    fn unit() -> Self;
    /// cast to f64
    fn to_f64(self) -> f64;
    /// cast to usize (by flooring)
    fn to_usize(self) -> usize;
    fn wrapping_add(self, rhs: Self) -> Self;
    fn wrapping_sub(self, rhs: Self) -> Self;
    fn large_const() -> Self;
    /// similary equal
    fn sim_eq(self, rhs: Self) -> bool;
    /// difference allowed to be regarded as a same value
    fn eps() -> Self;
}
impl FlowRateLike for usize {
    fn zero() -> usize {
        0
    }
    fn unit() -> usize {
        1
    }
    fn to_f64(self) -> f64 {
        self as f64
    }
    fn to_usize(self) -> usize {
        self
    }
    fn wrapping_add(self, rhs: Self) -> Self {
        self.wrapping_add(rhs)
    }
    fn wrapping_sub(self, rhs: Self) -> Self {
        self.wrapping_sub(rhs)
    }
    fn large_const() -> Self {
        100
    }
    fn sim_eq(self, rhs: Self) -> bool {
        // integer type does not need to consider the floating error
        self == rhs
    }
    fn eps() -> Self {
        0
    }
}
impl FlowRateLike for f64 {
    fn zero() -> Self {
        0.0
    }
    fn unit() -> Self {
        1.0
    }
    fn to_f64(self) -> f64 {
        self
    }
    fn to_usize(self) -> usize {
        // flooring
        self as usize
    }
    fn wrapping_add(self, rhs: Self) -> Self {
        // no overflow
        self + rhs
    }
    fn wrapping_sub(self, rhs: Self) -> Self {
        // no overflow
        self - rhs
    }
    fn large_const() -> Self {
        100.0
    }
    fn sim_eq(self, rhs: Self) -> bool {
        (self - rhs).abs() <= Self::eps()
    }
    fn eps() -> Self {
        0.000000001
    }
}

/// type of a cost (of edges per unit flow) in min-flow.
pub type Cost = f64;

//
// public functions
//

///
/// Find minimum cost flow on the FlowGraph
///
pub fn min_cost_flow<F, N, E>(graph: &DiGraph<N, E>) -> Option<Flow<F>>
where
    F: FlowRateLike,
    N: std::fmt::Debug,
    E: FlowEdge<F> + ConstCost + std::fmt::Debug,
{
    let init_flow = find_initial_flow(graph);

    match init_flow {
        Some(flow) => {
            draw_with_flow(graph, &flow);
            Some(min_cost_flow_from(graph, &flow))
        }
        None => None,
    }
}

///
/// Find minimum cost flow on the ConvexFlowGraph
///
pub fn min_cost_flow_convex<F, N, E>(graph: &DiGraph<N, E>) -> Option<Flow<F>>
where
    F: FlowRateLike,
    N: std::fmt::Debug,
    E: FlowEdge<F> + ConvexCost<F> + std::fmt::Debug,
{
    // (1) convert to normal FlowGraph and find the min-cost-flow
    let fg = match to_fixed_flow_graph(graph) {
        Some(fg) => fg,
        None => return None,
    };

    let fg_flow = match min_cost_flow(&fg) {
        Some(fg_flow) => fg_flow,
        None => return None,
    };

    // (2) convert-back to the flow on the ConvexFlowGraph
    Some(restore_convex_flow(&fg_flow, &fg, &graph))
}

///
/// Find minimum cost flow on the Graph whose edge is ConvexFlowEdge.
/// This solver requires less memory.
///
pub fn min_cost_flow_convex_fast<F, N, E>(graph: &DiGraph<N, E>) -> Option<Flow<F>>
where
    F: FlowRateLike,
    N: std::fmt::Debug,
    E: FlowEdge<F> + ConvexCost<F> + std::fmt::Debug,
{
    // (1) find the initial flow, by assigning constant cost to the flow.
    let init_flow = find_initial_flow(graph);

    // (2) upgrade the flow, by finding a negative cycle in residue graph.
    match init_flow {
        Some(flow) => {
            // draw_with_flow(graph, &flow);
            Some(min_cost_flow_from_convex(graph, &flow))
        }
        None => None,
    }
}

//
// internal functions
//

///
/// Find minimum cost flow of the special FlowGraph, whose demand is always zero.
///
pub fn min_cost_flow_from_zero<F, N, E>(graph: &DiGraph<N, E>) -> Flow<F>
where
    F: FlowRateLike,
    E: FlowEdge<F> + ConstCost,
{
    assert!(is_zero_demand_flow_graph(&graph));
    let flow = Flow::new(graph.edge_count(), F::zero());
    min_cost_flow_from(graph, &flow)
}

///
/// Find minimum cost by starting from the specified flow values.
///
pub fn min_cost_flow_from<F, N, E>(graph: &DiGraph<N, E>, init_flow: &Flow<F>) -> Flow<F>
where
    F: FlowRateLike,
    E: FlowEdge<F> + ConstCost,
{
    let mut flow = init_flow.clone();

    loop {
        assert_valid_flow(&flow, &graph);
        // solution of const cost is independent of cycle detection method
        // use BellmanFord because it is fastest.
        match improve_flow(graph, &flow, CycleDetectMethod::BellmanFord) {
            Some(new_flow) => {
                flow = new_flow;
                continue;
            }
            None => {
                break;
            }
        };
    }

    flow
}

///
/// Find minimum cost by starting from the specified flow values in ConvexCost Flowgraph.
///
pub fn min_cost_flow_from_convex<F, N, E>(graph: &DiGraph<N, E>, init_flow: &Flow<F>) -> Flow<F>
where
    F: FlowRateLike,
    E: FlowEdge<F> + ConvexCost<F>,
{
    let mut flow = init_flow.clone();

    // TODO assert graph edge has convex function?
    // assert!(is_convex_cost_flow_graph(graph));

    loop {
        assert!(is_valid_flow(&flow, &graph));
        // solution of convex cost is independent of cycle detection method
        // use BellmanFord because it is fastest.
        match improve_flow_convex(graph, &flow, CycleDetectMethod::BellmanFord) {
            Some(new_flow) => {
                flow = new_flow;
                continue;
            }
            None => {
                break;
            }
        };
    }

    flow
}

///
/// enumerate neighboring flows of current flow on MinFlowNetwork.
///
pub fn enumerate_neighboring_flows<F, N, E>(
    graph: &DiGraph<N, E>,
    flow: &Flow<F>,
    max_depth: Option<usize>,
) -> Vec<(Flow<F>, UpdateInfo)>
where
    F: FlowRateLike,
    E: FlowEdge<F> + ConvexCost<F>,
{
    let rg = flow_to_residue_convex(graph, flow);
    enumerate_neighboring_flows_in_residue(&rg, flow, max_depth)
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flow_rate_like_float() {
        assert_eq!(10.1_f64.to_usize(), 10);
        assert_eq!(10.9_f64.to_usize(), 10);
    }
}
