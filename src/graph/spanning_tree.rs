//!
//! Generate spanning tree by breadth first search from a specified node
//!
use super::cycle::Cycle;
use fnv::FnvHashSet as HashSet;
use itertools::Itertools;
use petgraph::algo::astar;
use petgraph::dot::Dot;
use petgraph::graph::{DefaultIx, DiGraph, EdgeIndex, Graph, NodeIndex, UnGraph};
use petgraph::visit::{Bfs, VisitMap, Visitable};
use std::collections::VecDeque;

#[derive(Clone, Debug)]
struct EdgeInfo {
    origin: EdgeIndex,
}

#[derive(Clone, Debug)]
struct SpanningTree {
    pub tree: UnGraph<(), EdgeInfo>,
    pub extra_edges: Vec<EdgeIndex>,
}

impl SpanningTree {
    fn n_cycle_basis(&self) -> usize {
        self.extra_edges.len()
    }
    ///
    /// get an `index`-th cycle basis (as an list of edges in the original graph)
    ///
    fn cycle_basis<N, E>(&self, graph: &UnGraph<N, E>, index: usize) -> Cycle {
        let extra_edge = self.extra_edges[index];
        let (v, w) = graph.edge_endpoints(extra_edge).unwrap();

        // shortest path between two nodes spanned by the extra_edge
        let (_, path) = astar(&self.tree, v, |n| n == w, |e| 1, |_| 0).expect("");

        // convert to edges in tree
        let edges_in_tree = nodes_to_edges(&self.tree, &path);

        // convert to edges in the original graph using the `origin` information
        let mut edges_in_graph: Vec<EdgeIndex> = edges_in_tree
            .into_iter()
            .map(|edge| self.tree.edge_weight(edge).unwrap().origin)
            .collect();
        edges_in_graph.push(extra_edge);

        Cycle::new(edges_in_graph)
    }
}

fn nodes_to_edges<N, E>(graph: &UnGraph<N, E>, nodes: &[NodeIndex]) -> Vec<EdgeIndex> {
    let mut edges = Vec::new();

    for (&v, &w) in nodes.iter().tuple_windows() {
        let edge = graph.find_edge(v, w).unwrap();
        edges.push(edge);
    }

    edges
}

///
/// create a spanning tree of Graph (direction of edges is ignored) by breadth first search from a
/// specified node.
///
/// ## Assumptions
///
/// - all nodes in the graph are reachable from the specified `start` node.
///
fn spanning_tree<N, E>(graph: &UnGraph<N, E>, start: NodeIndex) -> SpanningTree {
    let mut edges = Vec::new();
    let mut edge_used = vec![false; graph.edge_count()];

    // do bfs
    let mut visited = graph.visit_map();
    visited.visit(start);
    let mut queue = VecDeque::new();
    queue.push_back(start);

    while let Some(node) = queue.pop_front() {
        for succ in graph.neighbors(node) {
            if visited.visit(succ) {
                // if first visit
                queue.push_back(succ);

                // add this edge into spanning tree
                // parallel edges and self-loop is not allowed in tree
                // so only the first edge between two nodes is considered.
                let edge = graph.find_edge(node, succ).unwrap();
                edge_used[edge.index()] = true;
                edges.push((node, succ, EdgeInfo { origin: edge }));
            }
        }
    }

    let mut tree: UnGraph<(), EdgeInfo> = UnGraph::from_edges(&edges);
    println!("tree={:?}", Dot::with_config(&tree, &[]));

    // collect unused edges
    let mut extra_edges: Vec<EdgeIndex<DefaultIx>> = Vec::new();
    for (edge_index, is_used) in edge_used.into_iter().enumerate() {
        if !is_used {
            extra_edges.push(EdgeIndex::new(edge_index));
        }
    }

    println!("{:?}", extra_edges);

    SpanningTree { tree, extra_edges }
}

///
/// convert directed graph into undirected graph
///
fn to_undirected_graph<N, E>(digraph: DiGraph<N, E>) -> UnGraph<N, E> {
    digraph.into_edge_type()
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spanning_tree_mock_undirected_1() {
        let mut g: DiGraph<(), ()> = DiGraph::new();
        g.extend_with_edges(&[(0, 1), (2, 0), (2, 1)]);

        let mut g2 = to_undirected_graph(g);
        println!("{:?}", Dot::with_config(&g2, &[]));

        let st = spanning_tree(&g2, NodeIndex::new(0));
        for i in 0..st.n_cycle_basis() {
            println!("#{}", i);
            println!("{:?}", st.cycle_basis(&g2, i));
        }
    }

    #[test]
    fn spanning_tree_mock_undirected_2() {
        let mut g: UnGraph<(), ()> = UnGraph::from_edges(&[(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]);
        println!("{:?}", Dot::with_config(&g, &[]));

        let st = spanning_tree(&g, NodeIndex::new(0));
        for i in 0..st.n_cycle_basis() {
            println!("#{}", i);
            println!("{:?}", st.cycle_basis(&g, i));
        }
    }

    #[test]
    fn spanning_tree_mock_undirected_3() {
        let mut g: UnGraph<(), ()> = UnGraph::from_edges(&[
            (0, 1),
            (0, 2),
            (1, 2),
            (1, 3),
            (2, 3),
            (1, 3), // parallel edge
            (1, 1), // self loop
        ]);
        println!("{:?}", Dot::with_config(&g, &[]));

        let st = spanning_tree(&g, NodeIndex::new(0));
        for i in 0..st.n_cycle_basis() {
            println!("#{}", i);
            println!("{:?}", st.cycle_basis(&g, i));
        }
    }
}
