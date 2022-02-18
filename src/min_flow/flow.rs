//! Flow graph definitions
//! - FlowEdge, FlowEdgeRaw<T>
//! - FlowGraph, FlowGraphRaw<T>
//! - Flow
use crate::vector::{DenseStorage, EdgeVec};
use petgraph::graph::{DiGraph, EdgeIndex};
use petgraph::visit::EdgeRef; // for EdgeReference.id()
use petgraph::Direction;

/// Edge of FlowGraph
///
/// * `demand()`: demand `l(e)`
/// * `capacity()`: capacity `u(e)`
/// * `cost()`: cost per unit flow `c(e)`
///
/// `[l, u], c`
pub trait FlowEdgeLike {
    /// Demand of the edge, Lower limit of the flow
    fn demand(&self) -> u32;
    /// Capacity of the edge, Upper limit of the flow
    fn capacity(&self) -> u32;
    /// Cost-per-unit-flow of the edge
    fn cost(&self) -> f64;
}

/// Edge attributes used in FlowGraph
/// It has
/// - demand l
/// - capacity u
/// - cost per flow c
/// [l, u], c
///
/// it can contain additional information in T.
#[derive(Debug, Copy, Clone)]
pub struct FlowEdgeRaw<T> {
    /// demand (lower limit of flow) of the edge l(e)
    pub demand: u32,
    /// capacity (upper limit of flow) of the edge u(e)
    pub capacity: u32,
    /// cost per unit flow
    pub cost: f64,
    /// auxiliary informations
    pub info: T,
}

pub type FlowEdge = FlowEdgeRaw<()>;

impl FlowEdge {
    pub fn new(demand: u32, capacity: u32, cost: f64) -> FlowEdge {
        FlowEdge {
            demand,
            capacity,
            cost,
            info: (),
        }
    }
}

impl<T> std::fmt::Display for FlowEdgeRaw<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{},{}] {}", self.demand, self.capacity, self.cost)
    }
}

impl<T> FlowEdgeLike for FlowEdgeRaw<T> {
    fn demand(&self) -> u32 {
        self.demand
    }
    fn capacity(&self) -> u32 {
        self.capacity
    }
    fn cost(&self) -> f64 {
        self.cost
    }
}

/// FlowGraph definition
pub type FlowGraph = DiGraph<(), FlowEdge>;
pub type FlowGraphRaw<T> = DiGraph<(), FlowEdgeRaw<T>>;

/// Flow definitions
///
/// Flow f is a mapping of u32 f(e) to each edge e
pub type Flow = EdgeVec<DenseStorage<u32>>;

///
/// Check if the flow is valid, i.e. it satisfies
/// - flows of all edges are defined
/// - demand and capacity constraint
/// - flow constraint
///
pub fn is_valid_flow<N, E: FlowEdgeLike>(flow: &Flow, graph: &DiGraph<N, E>) -> bool {
    is_defined_for_all_edges(flow, graph)
        && is_in_demand_and_capacity(flow, graph)
        && is_satisfying_flow_constraint(flow, graph)
}

///
/// Check if the flow contains all edges
///
pub fn is_defined_for_all_edges<N, E: FlowEdgeLike>(flow: &Flow, graph: &DiGraph<N, E>) -> bool {
    flow.len() == graph.edge_count()
}

///
/// For each edge, the flow must satisfy `demand <= flow <= capacity`.
/// This function checks it
///
pub fn is_in_demand_and_capacity<N, E: FlowEdgeLike>(flow: &Flow, graph: &DiGraph<N, E>) -> bool {
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
pub fn is_satisfying_flow_constraint<N, E: FlowEdgeLike>(
    flow: &Flow,
    graph: &DiGraph<N, E>,
) -> bool {
    graph.node_indices().all(|v| {
        let in_flow: u32 = graph
            .edges_directed(v, Direction::Incoming)
            .map(|er| flow[er.id()])
            .sum();
        let out_flow: u32 = graph
            .edges_directed(v, Direction::Outgoing)
            .map(|er| flow[er.id()])
            .sum();
        in_flow == out_flow
    })
}

///
/// cost trait
///
pub trait EdgeCost {
    fn cost(&self, flow: u32) -> f64;
}

impl<T> EdgeCost for FlowEdgeRaw<T> {
    fn cost(&self, flow: u32) -> f64 {
        self.cost * flow as f64
    }
}

pub fn total_cost<N, E: EdgeCost>(graph: &DiGraph<N, E>, flow: &Flow) -> f64 {
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
