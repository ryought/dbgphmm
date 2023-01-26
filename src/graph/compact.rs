//!
//! Graph compacting functions
//! Functions in this file handles generic petgraph::DiGraph, not Dbg struct.
//!
//! * remove_deadends
//! * collapse_simple_paths
//!
use petgraph::{
    graph::{DiGraph, EdgeIndex, Graph, NodeIndex, UnGraph},
    stable_graph::StableDiGraph,
    visit::EdgeRef,
    Direction,
};
use std::collections::VecDeque;

///
/// get node degree (in degree if dir=Incoming, out degree if dir=Outgoing)
///
fn node_degree<N, E>(graph: &DiGraph<N, E>, node: NodeIndex, dir: Direction) -> usize {
    graph.edges_directed(node, dir).count()
}

///
/// get a single edge from/to a node.
///
fn find_edge_directed<N, E>(
    graph: &DiGraph<N, E>,
    node: NodeIndex,
    dir: Direction,
) -> Option<EdgeIndex> {
    graph.edges_directed(node, dir).next().map(|er| er.id())
}

///
/// find a deadend node in graph (in order of index)
///
fn find_deadend<N, E>(graph: &DiGraph<N, E>) -> Option<NodeIndex> {
    graph.node_indices().find(|&node| is_deadend(graph, node))
}

fn is_deadend<N, E>(graph: &DiGraph<N, E>, node: NodeIndex) -> bool {
    node_degree(graph, node, Direction::Incoming) == 0
        || node_degree(graph, node, Direction::Outgoing) == 0
}

fn is_deadend_stable<N, E>(graph: &StableDiGraph<N, E>, node: NodeIndex) -> bool {
    graph.edges_directed(node, Direction::Incoming).count() == 0
        || graph.edges_directed(node, Direction::Outgoing).count() == 0
}

///
/// Remove all deadend nodes
/// * in_degree is zero
/// * out_degree is zero
///
/// ## Notes
/// * Removing can cause node index change
/// * There is room to improve, current algorithm takes O(N^2).
///
pub fn remove_deadends<N, E>(graph: &mut DiGraph<N, E>) {
    while let Some(deadend) = find_deadend(&graph) {
        graph.remove_node(deadend);
    }
}

///
/// Remove all deadend nodes (in_deg == out_deg == 0) faster than O(V^2)
///
/// ## Notes
/// * Removing can cause node index change
///
pub fn remove_deadends_fast<N, E>(graph: DiGraph<N, E>) -> DiGraph<N, E> {
    //
    let mut g: StableDiGraph<N, E> = graph.into();

    let mut queue: VecDeque<NodeIndex> = g
        .node_indices()
        .filter(|&node| is_deadend_stable(&g, node))
        .collect();

    while let Some(deadend) = queue.pop_front() {
        if g.contains_node(deadend) && is_deadend_stable(&g, deadend) {
            // child of deadend is added to the queue
            for child in g.neighbors_directed(deadend, Direction::Incoming) {
                queue.push_back(child);
            }
            for child in g.neighbors_directed(deadend, Direction::Outgoing) {
                queue.push_back(child);
            }
            g.remove_node(deadend);
        }
    }

    // convert to normal digraph
    g.into()
}

///
/// find a node in simple-path i.e. in-deg = out-deg = 1
///
/// Except nodes with a single self-loop.
///
fn find_simple_path_node<N, E>(graph: &DiGraph<N, E>) -> Option<NodeIndex> {
    graph.node_indices().find(|&node| {
        node_degree(graph, node, Direction::Incoming) == 1
            && node_degree(graph, node, Direction::Outgoing) == 1
            && find_edge_directed(graph, node, Direction::Incoming)
                != find_edge_directed(graph, node, Direction::Outgoing)
    })
}

///
/// remove all nodes whose in/out-degree is 1
///
pub fn compact_simple_paths<N: Clone, E: Clone>(
    graph: &DiGraph<N, E>,
) -> DiGraph<N, Vec<(EdgeIndex, E)>> {
    let mut ret = graph.map(
        |_node, node_weight| node_weight.clone(),
        |edge, edge_weight| vec![(edge, edge_weight.clone())],
    );

    while let Some(node) = find_simple_path_node(&ret) {
        let in_edge_ref = ret
            .edges_directed(node, Direction::Incoming)
            .next()
            .unwrap();
        let out_edge_ref = ret
            .edges_directed(node, Direction::Outgoing)
            .next()
            .unwrap();
        let s = in_edge_ref.source();
        let t = out_edge_ref.target();
        let in_edge = in_edge_ref.id();
        let out_edge = out_edge_ref.id();
        if in_edge == out_edge {
            // self loop
        }

        // take edge weight
        // removing in_edge changes out_edge's index
        // node with larger index should be removed first.
        let (in_edge_weight, out_edge_weight) = if in_edge.index() < out_edge.index() {
            // out fisrt
            let o = ret.remove_edge(out_edge).unwrap();
            let i = ret.remove_edge(in_edge).unwrap();
            (i, o)
        } else {
            // in first
            let i = ret.remove_edge(in_edge).unwrap();
            let o = ret.remove_edge(out_edge).unwrap();
            (i, o)
        };
        let edge_weight = [in_edge_weight, out_edge_weight].concat();

        // removing node causes index change
        // so add_edge must be done first
        ret.add_edge(s, t, edge_weight);
        ret.remove_node(node);
    }

    ret
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{ei, ni};
    use crate::graph::utils::{to_edge_list, to_node_list};
    use petgraph::dot::Dot;

    #[test]
    fn compact_simple_paths_test_1() {
        let g: DiGraph<(), usize> = DiGraph::from_edges(&[
            (0, 1, 10),
            (1, 2, 20),
            (2, 4, 30),
            (0, 3, 40),
            (3, 4, 50),
            (4, 0, 60),
        ]);
        println!("{:?}", Dot::with_config(&g, &[]));
        let h = compact_simple_paths(&g);
        println!("{:?}", Dot::with_config(&h, &[]));
        assert_eq!(to_node_list(&h), vec![(0, ()), (1, ())]);
        assert_eq!(
            to_edge_list(&h),
            vec![
                (0, 0, 1, vec![(ei(0), 10), (ei(1), 20), (ei(2), 30)]),
                (1, 1, 0, vec![(ei(5), 60)]),
                (2, 0, 1, vec![(ei(3), 40), (ei(4), 50)]),
            ]
        );
    }

    #[test]
    fn remove_deadends_test() {
        // graph with deadends
        {
            let mut g: DiGraph<(), ()> = DiGraph::from_edges(&[
                // main circle <1,2,3,4>
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 1),
                //deadend edges 0,5,6
                (2, 0),
                (4, 5),
                (6, 5),
            ]);
            println!("{:?}", Dot::with_config(&g, &[]));
            assert_eq!(node_degree(&g, ni(4), Direction::Incoming), 1);
            assert_eq!(node_degree(&g, ni(4), Direction::Outgoing), 2);
            assert_eq!(node_degree(&g, ni(0), Direction::Incoming), 1);
            assert_eq!(node_degree(&g, ni(0), Direction::Outgoing), 0);
            assert_eq!(find_deadend(&g), Some(ni(0)));

            remove_deadends(&mut g);
            println!("{:?}", Dot::with_config(&g, &[]));
            assert_eq!(g.node_count(), 4);
            assert_eq!(g.edge_count(), 4);

            let g2 = remove_deadends_fast(g.clone());
            println!("{:?}", Dot::with_config(&g2, &[]));
            assert_eq!(g2.node_count(), 4);
            assert_eq!(g2.edge_count(), 4);
        }

        {
            let mut g: DiGraph<(), ()> = DiGraph::from_edges(&[
                // main circle <1,2,3,4>
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 1),
                // deadend nodes 0,8 5,6,7
                (2, 0),
                (0, 8),
                (4, 5),
                (6, 5),
                (5, 7),
            ]);
            println!("{:?}", Dot::with_config(&g, &[]));

            remove_deadends(&mut g);
            println!("{:?}", Dot::with_config(&g, &[]));
            assert_eq!(g.node_count(), 4);
            assert_eq!(g.edge_count(), 4);

            let g2 = remove_deadends_fast(g.clone());
            println!("{:?}", Dot::with_config(&g2, &[]));
            assert_eq!(g2.node_count(), 4);
            assert_eq!(g2.edge_count(), 4);
        }

        // graph with deadends
        {
            let mut g: DiGraph<(), ()> = DiGraph::from_edges(&[(0, 1), (1, 2), (2, 3), (3, 0)]);
            println!("{:?}", Dot::with_config(&g, &[]));
            assert_eq!(find_deadend(&g), None);

            remove_deadends(&mut g);
            println!("{:?}", Dot::with_config(&g, &[]));
            assert_eq!(g.node_count(), 4);
            assert_eq!(g.edge_count(), 4);

            let g2 = remove_deadends_fast(g.clone());
            println!("{:?}", Dot::with_config(&g2, &[]));
            assert_eq!(g2.node_count(), 4);
            assert_eq!(g2.edge_count(), 4);
        }
    }
}
