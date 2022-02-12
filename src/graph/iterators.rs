//!
//! Extension to `petgraph` basic iterators
//!
use petgraph::graph::DiGraph;
use petgraph::graph::{EdgeIndex, NodeIndex};
use petgraph::graph::{EdgeReferences, Edges, NodeReferences};
use petgraph::visit::EdgeRef;
use petgraph::visit::IntoNodeReferences;
use petgraph::{Directed, Direction};

/// Iterator struct for `nodes()`
/// (with `NodeReferences` of petgraph)
///
pub struct NodesIterator<'a, N: 'a> {
    nodes: NodeReferences<'a, N>,
}

impl<'a, N> NodesIterator<'a, N> {
    /// Create NodesIterator from the reference of DiGraph
    pub fn new<E>(graph: &'a DiGraph<N, E>) -> Self {
        NodesIterator {
            nodes: graph.node_references(),
        }
    }
}

impl<'a, N> Iterator for NodesIterator<'a, N> {
    type Item = (NodeIndex, &'a N);
    fn next(&mut self) -> Option<Self::Item> {
        self.nodes.next()
    }
}

/// Iterator struct for `edges()`
/// (with `NodeReferences` of petgraph)
///
pub struct EdgesIterator<'a, E: 'a> {
    edges: EdgeReferences<'a, E>,
}

impl<'a, E> EdgesIterator<'a, E> {
    /// Create EdgesIterator from the reference of DiGraph
    pub fn new<N>(graph: &'a DiGraph<N, E>) -> Self {
        EdgesIterator {
            edges: graph.edge_references(),
        }
    }
}

impl<'a, E> Iterator for EdgesIterator<'a, E> {
    type Item = (EdgeIndex, NodeIndex, NodeIndex, &'a E);
    fn next(&mut self) -> Option<Self::Item> {
        // extract edge reference
        match self.edges.next() {
            Some(er) => Some((er.id(), er.source(), er.target(), er.weight())),
            None => None,
        }
    }
}

/// Iterator for `childs()`
pub struct ChildEdges<'a, E: 'a> {
    edges: Edges<'a, E, Directed>,
}

impl<'a, E> ChildEdges<'a, E> {
    /// Create ChildEdges from the reference of DiGraph
    pub fn new<N>(graph: &'a DiGraph<N, E>, node: NodeIndex) -> Self {
        ChildEdges {
            edges: graph.edges_directed(node, Direction::Outgoing),
        }
    }
}

impl<'a, E> Iterator for ChildEdges<'a, E> {
    type Item = (EdgeIndex, NodeIndex, &'a E);
    fn next(&mut self) -> Option<Self::Item> {
        // edge reference
        match self.edges.next() {
            // er.source() = the given node
            // er.target() = child
            Some(er) => Some((er.id(), er.target(), er.weight())),
            None => None,
        }
    }
}

/// Iterator for `parents()`
pub struct ParentEdges<'a, E: 'a> {
    pub edges: Edges<'a, E, Directed>,
}

impl<'a, E> ParentEdges<'a, E> {
    /// Create ParentEdges from the reference of DiGraph
    pub fn new<N>(graph: &'a DiGraph<N, E>, node: NodeIndex) -> Self {
        ParentEdges {
            edges: graph.edges_directed(node, Direction::Incoming),
        }
    }
}

impl<'a, E> Iterator for ParentEdges<'a, E> {
    type Item = (EdgeIndex, NodeIndex, &'a E);
    fn next(&mut self) -> Option<Self::Item> {
        // edge reference
        match self.edges.next() {
            // er.source() = parent
            // er.target() = the given node
            Some(er) => Some((er.id(), er.source(), er.weight())),
            None => None,
        }
    }
}
