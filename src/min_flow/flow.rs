//! Flow graph definitions
//! - FlowEdge, FlowEdgeRaw<T>
//! - FlowGraph, FlowGraphRaw<T>
//! - Flow
use super::convex::ConvexCost;
use super::{ConstCost, Cost, FlowEdge, FlowRateLike};
use crate::vector::{DenseStorage, EdgeVec};
use petgraph::graph::DiGraph;
use petgraph::visit::EdgeRef; // for EdgeReference.id()
use petgraph::Direction;

/// Flow definitions
///
/// Flow f is a mapping of FlowRate(u32) f(e) to each edge e
pub type Flow<F> = EdgeVec<DenseStorage<F>>;

///
/// Check if the flow is valid, i.e. it satisfies
/// - flows of all edges are defined
/// - demand and capacity constraint
/// - flow constraint
///
pub fn is_valid_flow<F: FlowRateLike, N, E: FlowEdge<F>>(
    flow: &Flow<F>,
    graph: &DiGraph<N, E>,
) -> bool {
    is_defined_for_all_edges(flow, graph)
        && is_in_demand_and_capacity(flow, graph)
        && is_satisfying_flow_constraint(flow, graph)
}

pub fn assert_valid_flow<F: FlowRateLike, N, E: FlowEdge<F>>(
    flow: &Flow<F>,
    graph: &DiGraph<N, E>,
) {
    assert!(
        is_defined_for_all_edges(flow, graph),
        "flow does not match the graph size"
    );
    assert!(
        is_in_demand_and_capacity(flow, graph),
        "flow not in [demand,capacity]"
    );
    assert!(
        is_satisfying_flow_constraint(flow, graph),
        "flow does not satisfy flow constraint"
    );
}

///
/// Check if the flow contains all edges
///
pub fn is_defined_for_all_edges<F: FlowRateLike, N, E: FlowEdge<F>>(
    flow: &Flow<F>,
    graph: &DiGraph<N, E>,
) -> bool {
    flow.len() == graph.edge_count()
}

///
/// For each edge, the flow must satisfy `demand <= flow <= capacity`.
/// This function checks it
///
pub fn is_in_demand_and_capacity<F: FlowRateLike, N, E: FlowEdge<F>>(
    flow: &Flow<F>,
    graph: &DiGraph<N, E>,
) -> bool {
    graph.edge_indices().all(|e| {
        let ew = graph.edge_weight(e).unwrap();
        let f = flow[e];
        (ew.demand() <= f) && (f <= ew.capacity())
    })
}

///
/// For each node,
/// (the sum of out-going flows) should be equal to (the sum of in-coming flows).
///
pub fn is_satisfying_flow_constraint<F: FlowRateLike, N, E: FlowEdge<F>>(
    flow: &Flow<F>,
    graph: &DiGraph<N, E>,
) -> bool {
    graph.node_indices().all(|v| {
        let in_flow: F = graph
            .edges_directed(v, Direction::Incoming)
            .map(|er| flow[er.id()])
            .sum();
        let out_flow: F = graph
            .edges_directed(v, Direction::Outgoing)
            .map(|er| flow[er.id()])
            .sum();
        in_flow.sim_eq(out_flow)
    })
}

///
/// inspect flows by each node for debugging
///
pub fn inspect_flow_constraint<F: FlowRateLike, N, E: FlowEdge<F>>(
    flow: &Flow<F>,
    graph: &DiGraph<N, E>,
) {
    println!("[flow_inspector] start");
    for node in graph.node_indices() {
        println!("node#{}", node.index());
        for e in graph.edges_directed(node, Direction::Incoming) {
            println!("edge+ #{} {}", e.id().index(), flow[e.id()]);
        }
        for e in graph.edges_directed(node, Direction::Outgoing) {
            println!("edge- #{} {}", e.id().index(), flow[e.id()]);
        }
    }
    println!("[flow_inspector] end");
}

///
/// cost trait
///
pub trait EdgeCost<F: FlowRateLike> {
    fn cost(&self, flow: F) -> Cost;
}

impl<F: FlowRateLike, E: ConstCost + FlowEdge<F>> ConvexCost<F> for E {
    fn convex_cost(&self, flow: F) -> Cost {
        self.cost() * flow.to_f64()
    }
}

impl<F: FlowRateLike, E: ConvexCost<F>> EdgeCost<F> for E {
    fn cost(&self, flow: F) -> Cost {
        self.convex_cost(flow)
    }
}

///
/// Calculate the total cost of the flow in the graph.
///
pub fn total_cost<F: FlowRateLike, N, E: EdgeCost<F>>(
    graph: &DiGraph<N, E>,
    flow: &Flow<F>,
) -> Cost {
    graph
        .edge_indices()
        .map(|e| {
            let ew = graph.edge_weight(e).unwrap();
            let f = flow[e];
            ew.cost(f)
        })
        .sum()
}

//
// tests
//
#[cfg(test)]
mod tests {
    use super::super::mocks::mock_flow_network1;
    use super::super::utils::draw;
    use super::*;
    use petgraph::graph::EdgeIndex;

    #[test]
    fn flow_valid_tests() {
        let (g, _) = mock_flow_network1();
        draw(&g);

        // this is valid flow
        let f1 = Flow::from_vec(
            3,
            0,
            &[
                (EdgeIndex::new(0), 5),
                (EdgeIndex::new(1), 5),
                (EdgeIndex::new(2), 5),
            ],
        );
        assert!(is_defined_for_all_edges(&f1, &g));
        assert!(is_in_demand_and_capacity(&f1, &g));
        assert!(is_satisfying_flow_constraint(&f1, &g));
        assert!(is_valid_flow(&f1, &g));

        // this flow overs the capacity
        let f2 = Flow::from_vec(
            3,
            0,
            &[
                (EdgeIndex::new(0), 100),
                (EdgeIndex::new(1), 100),
                (EdgeIndex::new(2), 100),
            ],
        );
        assert!(is_defined_for_all_edges(&f2, &g));
        assert!(!is_in_demand_and_capacity(&f2, &g));
        assert!(is_satisfying_flow_constraint(&f2, &g));
        assert!(!is_valid_flow(&f2, &g));

        // this is a flow which not satisfies the flow constraint
        let f3 = Flow::from_vec(
            3,
            0,
            &[
                (EdgeIndex::new(0), 1),
                (EdgeIndex::new(1), 5),
                (EdgeIndex::new(2), 1),
            ],
        );
        assert!(is_defined_for_all_edges(&f3, &g));
        assert!(is_in_demand_and_capacity(&f3, &g));
        assert!(!is_satisfying_flow_constraint(&f3, &g));
        assert!(!is_valid_flow(&f3, &g));

        // this is a partial flow
        let f4 = Flow::from_vec(1, 0, &[(EdgeIndex::new(0), 1)]);
        assert!(!is_defined_for_all_edges(&f4, &g));
        assert!(!is_valid_flow(&f4, &g));
    }
}
