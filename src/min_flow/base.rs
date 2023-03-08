//!
//! # `FlowGraph`
//!
//! Example implementation of flow network
//!
use super::{ConstCost, Cost, FlowEdge, FlowRateLike};
use petgraph::graph::DiGraph;

/// FlowGraph definition
pub type FlowGraph<F> = DiGraph<(), FlowEdgeBase<F>>;
pub type FlowGraphRaw<F, T> = DiGraph<(), FlowEdgeRaw<F, T>>;

/// Edge attributes used in FlowGraph
/// It has
/// - demand l
/// - capacity u
/// - cost per flow c
/// [l, u], c
///
/// it can contain additional information in T.
#[derive(Debug, Copy, Clone)]
pub struct FlowEdgeRaw<F: FlowRateLike, T> {
    /// demand (lower limit of flow) of the edge l(e)
    pub demand: F,
    /// capacity (upper limit of flow) of the edge u(e)
    pub capacity: F,
    /// cost per unit flow
    pub cost: Cost,
    /// auxiliary informations
    pub info: T,
}

pub type FlowEdgeBase<F> = FlowEdgeRaw<F, ()>;

impl<F: FlowRateLike> FlowEdgeBase<F> {
    pub fn new(demand: F, capacity: F, cost: Cost) -> FlowEdgeBase<F> {
        FlowEdgeBase {
            demand,
            capacity,
            cost,
            info: (),
        }
    }
}

impl<F: FlowRateLike + std::fmt::Display, T> std::fmt::Display for FlowEdgeRaw<F, T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{},{}] {}", self.demand, self.capacity, self.cost)
    }
}

impl<F: FlowRateLike, T> FlowEdge<F> for FlowEdgeRaw<F, T> {
    fn demand(&self) -> F {
        self.demand
    }
    fn capacity(&self) -> F {
        self.capacity
    }
}

impl<F: FlowRateLike, T> ConstCost for FlowEdgeRaw<F, T> {
    fn cost(&self) -> Cost {
        self.cost
    }
}
