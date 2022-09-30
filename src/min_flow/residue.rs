//! Residue graph related definitions
//! - ResidueEdge
//! - ResidueGraph
//! - ResidueDirection
//!
use super::convex::ConvexCost;
use super::flow::{ConstCost, EdgeCost, Flow, FlowEdge};
use super::utils::draw;
use super::{Cost, FlowRate};
use crate::graph::bellman_ford::HasEpsilon;
use crate::graph::float_weight::{
    edge_cycle_to_node_cycle, is_cycle, is_edge_simple, is_negative_cycle, node_list_to_edge_list,
    total_weight,
};
use crate::graph::min_mean_weight_cycle::edge_cond::find_negative_cycle_with_edge_cond;
use crate::graph::min_mean_weight_cycle::{find_negative_cycle, find_negative_edge_cycle};
use crate::graph::FloatWeight;
use itertools::Itertools; // for tuple_windows
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
use petgraph::prelude::*;
use petgraph::visit::VisitMap;
use std::cmp::Ordering;

// basic definitions

/// Edge attributes used in ResidueGraph
#[derive(Debug, Default, Copy, Clone)]
pub struct ResidueEdge {
    /// The movable amount of the flow
    pub count: FlowRate,
    /// Cost of the unit change of this flow
    pub weight: Cost,
    /// Original edge index of the source graph
    pub target: EdgeIndex,
    /// +1 or -1
    pub direction: ResidueDirection,
}

impl ResidueEdge {
    pub fn new(
        count: FlowRate,
        weight: Cost,
        target: EdgeIndex,
        direction: ResidueDirection,
    ) -> ResidueEdge {
        ResidueEdge {
            count,
            weight,
            target,
            direction,
        }
    }
    pub fn only_weight(weight: Cost) -> ResidueEdge {
        ResidueEdge {
            weight,
            // filled by default values
            target: EdgeIndex::new(0),
            count: 0,
            direction: ResidueDirection::Up,
        }
    }
}

impl FloatWeight for ResidueEdge {
    fn float_weight(&self) -> f64 {
        self.weight
    }
    fn epsilon() -> f64 {
        0.00001
    }
}

/// Residue direction enum
/// residue edge has two types
#[derive(Debug, Copy, Clone, PartialEq)]
pub enum ResidueDirection {
    /// Up edge: it can increase(+1) flow
    Up,
    /// Down edge: it can decrease(-1) flow
    Down,
}

impl Default for ResidueDirection {
    fn default() -> Self {
        ResidueDirection::Up
    }
}

/// ResidueGraph definition
pub type ResidueGraph = DiGraph<(), ResidueEdge>;

//
// conversion functions
//

/// Convert FlowGraph with Flow into ResidueGraph.
///
/// FlowGraph and Flow
/// v -> w
///  e = ([l,u],c), f
///
/// into
///
/// ResidueGraph
/// v -> w
///  e1 = (u-f, +c) if u-f>0
/// w -> v
///  e2 = (f-l, -c) if f-l>0
pub fn flow_to_residue<N, E: FlowEdge + ConstCost>(
    graph: &DiGraph<N, E>,
    flow: &Flow,
) -> ResidueGraph {
    let mut rg: ResidueGraph = ResidueGraph::new();

    // create two edges (Up and Down) for each edge
    for e in graph.edge_indices() {
        let f = flow[e];
        let ew = graph.edge_weight(e).unwrap();
        let (v, w) = graph.edge_endpoints(e).unwrap();

        let mut edges = Vec::new();
        if f < ew.capacity() {
            // up movable
            edges.push((
                v,
                w,
                ResidueEdge::new(ew.capacity() - f, ew.cost(), e, ResidueDirection::Up),
            ));
        }
        if f > ew.demand() {
            // down movable
            edges.push((
                w,
                v,
                ResidueEdge::new(f - ew.demand(), -ew.cost(), e, ResidueDirection::Down),
            ));
        }
        rg.extend_with_edges(&edges);
    }
    rg
}

/// Convert FlowGraph with Flow with ConvexCost into ResidueGraph.
///
/// For each edge in FlowGraph with Flow
/// ```text
/// e(v -> w) = ([l,u],c), f
/// ```
///
/// create two edges in ResidueGraph
/// ```text
/// e1(v -> w) = (1, c(f+1) - c(f)) if u - f > 0
///
/// e2(w -> v) = (1, c(f-1) - c(f)) if f - l > 0
/// ```
pub fn flow_to_residue_convex<N, E>(graph: &DiGraph<N, E>, flow: &Flow) -> ResidueGraph
where
    E: FlowEdge + ConvexCost,
{
    let mut rg: ResidueGraph = ResidueGraph::new();

    // create two edges (Up and Down) for each edge
    for e in graph.edge_indices() {
        let f = flow[e];
        let ew = graph.edge_weight(e).unwrap();
        let (v, w) = graph.edge_endpoints(e).unwrap();

        let mut edges = Vec::new();
        if f < ew.capacity() {
            // up movable
            // cost=cost(f+1)-cost(f)
            edges.push((
                v,
                w,
                ResidueEdge::new(1, ew.cost_diff(f), e, ResidueDirection::Up),
            ));
        }
        if f > ew.demand() {
            // down movable
            // cost=cost(f-1)-cost(f)
            edges.push((
                w,
                v,
                ResidueEdge::new(1, -ew.cost_diff(f - 1), e, ResidueDirection::Down),
            ));
        }

        // if up/down movable,
        // self round loop (v->w and w->v) should not have negative weight.
        if f < ew.capacity() && f > ew.demand() {
            let cost_up = ew.cost_diff(f);
            let cost_down = -ew.cost_diff(f - 1);

            // TODO this assertion is valid only if the cost function is convex.
            // assert!(cost_up + cost_down >= 0.0);
        }

        rg.extend_with_edges(&edges);
    }
    rg
}

#[allow(dead_code)]
fn residue_to_float_weighted_graph(graph: &ResidueGraph) -> DiGraph<(), Cost> {
    graph.map(|_, _| (), |_, ew| ew.weight)
}

//
// internal functions to find a update of the flow
// (i.e. the negative cycle in ResidueGraph)
//

///
/// Update the flow by a negative cycle on a residue graph.
///
fn apply_residual_edges_to_flow(flow: &Flow, rg: &ResidueGraph, edges: &[EdgeIndex]) -> Flow {
    let mut new_flow = flow.clone();

    // (1) determine flow_change_amount
    // that is the minimum of ResidueEdge.count
    let flow_change_amount = edges
        .iter()
        .map(|&e| {
            let ew = rg.edge_weight(e).unwrap();
            ew.count
        })
        .min()
        .unwrap();

    // (2) apply these changes to the flow along the cycle
    for edge in edges {
        let ew = rg.edge_weight(*edge).unwrap();
        // convert back to the original edgeindex
        let original_edge = ew.target;

        // use `wrapping_{add,sub}` because
        // in the some ordering of residue edges, applying -1 on a zero-flow edge can happen.
        // As long as the residue edges is valid (i.e. it makes cycle in the residue graph)
        // the final flow should satisfy the flow condition.
        new_flow[original_edge] = match ew.direction {
            ResidueDirection::Up => new_flow[original_edge].wrapping_add(flow_change_amount),
            ResidueDirection::Down => new_flow[original_edge].wrapping_sub(flow_change_amount),
        };
    }

    new_flow
}

fn find_negative_cycle_in_whole_graph(graph: &ResidueGraph) -> Option<Vec<EdgeIndex>> {
    let mut node = NodeIndex::new(0);
    let mut dfs = Dfs::new(&graph, node);

    loop {
        // find negative cycle with prohibiting e1 -> e2 transition
        //
        //    e1 (+1 of e)
        // v --->
        //   <--- w
        //    e2 (-1 of e)
        //
        let path = find_negative_cycle_with_edge_cond(&graph, node, |e_a, e_b| {
            let ew_a = graph.edge_weight(e_a).unwrap();
            let ew_b = graph.edge_weight(e_b).unwrap();

            let target_is_different = ew_a.target != ew_b.target;
            let dir_is_same = ew_a.direction == ew_b.direction;
            target_is_different || dir_is_same
        });

        if path.is_some() {
            return path;
        }

        // search for alternative start point
        dfs.move_to(node);
        while let Some(_nx) = dfs.next(&graph) {}
        let unvisited_node = graph
            .node_indices()
            .find(|node| !dfs.discovered.is_visited(node));

        // if there is unvisited node, search again for negative cycle
        match unvisited_node {
            Some(n) => {
                node = n;
                continue;
            }
            None => break,
        }
    }
    return None;
}

///
/// Update residue graph by finding negative cycle
///
pub fn improve_residue_graph(rg: &ResidueGraph) -> Option<Vec<EdgeIndex>> {
    // find negative weight cycles
    let path = find_negative_cycle_in_whole_graph(&rg);
    // draw(&rg);

    match path {
        Some(edges) => {
            // check if this is actually negative cycle
            assert!(
                is_cycle(&rg, &edges),
                "the cycle was not valid. edges={:?}",
                edges
            );
            assert!(
                is_edge_simple(&rg, &edges),
                "the cycle is not edge-simple. edges={:?}",
                edges
            );
            assert!(
                is_negative_cycle(&rg, &edges),
                "total weight of the found negative cycle is not negative. edges={:?} total_weight={}",
                edges,
                total_weight(&rg, &edges)
            );

            // TODO assert using is_meaningful_cycle?

            Some(edges)
        }
        None => None,
    }
}

///
/// WIP
///
fn is_meaningful_cycle(rg: &ResidueGraph, cycle: &[EdgeIndex]) -> bool {
    unimplemented!();
}

/// create a new improved flow from current flow
/// by upgrading along the negative weight cycle in the residual graph
fn update_flow_in_residue_graph(flow: &Flow, rg: &ResidueGraph) -> Option<(Flow, Vec<EdgeIndex>)> {
    match improve_residue_graph(rg) {
        Some(cycle) => {
            // apply these changes along the cycle to current flow
            let new_flow = apply_residual_edges_to_flow(&flow, &rg, &cycle);

            // if applying edges did not changed the flow (i.e. the edges was meaningless)
            // improve should fail.
            if &new_flow == flow {
                println!("meaningless cycle was found!");
                None
            } else {
                Some((new_flow, cycle))
            }
        }
        None => None,
    }
}

fn cycle_in_residue_graph_into_update_info(rg: &ResidueGraph, cycle: &[EdgeIndex]) -> UpdateInfo {
    cycle
        .iter()
        .map(|&e| {
            let ew = rg.edge_weight(e).unwrap();
            (ew.target, ew.direction)
        })
        .collect()
}

//
// public functions
//

/// create a new improved flow from current flow
/// by upgrading along the negative weight cycle in the residual graph
pub fn improve_flow<N, E: FlowEdge + ConstCost>(
    graph: &DiGraph<N, E>,
    flow: &Flow,
) -> Option<Flow> {
    let rg = flow_to_residue(graph, flow);
    match update_flow_in_residue_graph(flow, &rg) {
        Some((new_flow, _)) => Some(new_flow),
        None => None,
    }
}

///
/// information of updating a edge of either direction?
///
pub type UpdateInfo = Vec<(EdgeIndex, ResidueDirection)>;

///
/// summary of UpdateInfo
///
pub type UpdateSummary = Vec<(Vec<EdgeIndex>, ResidueDirection)>;

///
/// Convert a update cycle
///     [(1, +), (2, +), (3, +), (3, -), (2, -), (1, -)]
///     [(2, +), (3, +), (3, -), (2, -), (1, -), (1, +)]
/// into a normalized summary
///     [([1,2,3], +), ([3,2,1], -)]
///
fn to_contiguous_direction_list(
    updates: &[(EdgeIndex, ResidueDirection)],
) -> Vec<(Vec<EdgeIndex>, ResidueDirection)> {
    unimplemented!();
}

/// create a new improved flow from current flow
/// by upgrading along the negative weight cycle in the residual graph
pub fn improve_flow_convex_with_update_info<N, E>(
    graph: &DiGraph<N, E>,
    flow: &Flow,
) -> Option<(Flow, UpdateInfo)>
where
    E: FlowEdge + ConvexCost,
{
    let rg = flow_to_residue_convex(graph, flow);
    match update_flow_in_residue_graph(flow, &rg) {
        Some((new_flow, cycle)) => Some((
            new_flow,
            cycle_in_residue_graph_into_update_info(&rg, &cycle),
        )),
        None => None,
    }
}

/// create a new improved flow from current flow
/// by upgrading along the negative weight cycle in the residual graph
pub fn improve_flow_convex<N, E>(graph: &DiGraph<N, E>, flow: &Flow) -> Option<Flow>
where
    E: FlowEdge + ConvexCost,
{
    let rg = flow_to_residue_convex(graph, flow);
    match update_flow_in_residue_graph(flow, &rg) {
        Some((new_flow, _)) => Some(new_flow),
        None => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ei;

    #[test]
    fn petgraph_negative_cycle_test() {
        // small cycle test
        let mut g: DiGraph<(), f64> = Graph::new();
        let a = g.add_node(());
        let b = g.add_node(());
        g.add_edge(a, b, -10.0);
        g.add_edge(b, a, 9.0);
        let path = find_negative_cycle(&g, NodeIndex::new(0));
        assert_eq!(path.is_some(), true);
        let nodes = path.unwrap();
        assert!(nodes.contains(&NodeIndex::new(0)));
        assert!(nodes.contains(&NodeIndex::new(1)));
    }

    #[test]
    fn petgraph_negative_cycle_test2() {
        // self loop test, it will work fine
        let mut g: DiGraph<(), f64> = Graph::new();
        let a = g.add_node(());
        g.add_edge(a, a, -10.0);
        let path = find_negative_cycle(&g, NodeIndex::new(0));
        assert_eq!(path.is_some(), true);
        let nodes = path.unwrap();
        assert!(nodes.contains(&NodeIndex::new(0)));
    }

    #[test]
    fn negative_cycle_in_whole() {
        let mut g: ResidueGraph = ResidueGraph::new();
        let a = g.add_node(());
        let b = g.add_node(());
        let c = g.add_node(());
        g.add_edge(
            a,
            b,
            ResidueEdge::new(1, 10.0, EdgeIndex::new(0), ResidueDirection::Up),
        );
        g.add_edge(
            b,
            a,
            ResidueEdge::new(1, -1.0, EdgeIndex::new(1), ResidueDirection::Up),
        );
        g.add_edge(
            c,
            c,
            ResidueEdge::new(1, -1.0, EdgeIndex::new(2), ResidueDirection::Up),
        );
        let path = find_negative_cycle_in_whole_graph(&g);
        assert_eq!(path.is_some(), true);
        assert_eq!(path, Some(vec![ei(2)]));
    }
}
