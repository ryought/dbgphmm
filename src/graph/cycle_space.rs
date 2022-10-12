//!
//! cycle space
//!
//! enumerate all cycles in the undirected graph
//!
use super::cycle::{Cycle, SimpleCycle};
use fnv::FnvHashMap as HashMap;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex, UnGraph};
use petgraph::visit::EdgeRef;
use petgraph::Direction;
// unused imports
// use petgraph::unionfind::UnionFind;

#[derive(Clone, Debug)]
struct CycleSpace {
    ///
    /// a set of all cycle basis.
    ///
    basis: Vec<SimpleCycle>,
}

///
/// convert a graph `G` into adjacent graph `G*`
///
/// - an edge in `G` corresponds to a node in `G*`
/// - a pair of nodes (v, w) in `G*` will be connected by an edge v->w if corresponding edges in `G`
///   if the corresponding edge ends/starts from the same node.
///
fn into_adj_graph<N, F>(graph: &DiGraph<N, ()>, node_map: F) -> DiGraph<N, ()>
where
    F: Fn(&N, &N) -> N,
{
    let mut g = DiGraph::new();
    // map from edge in G into node in G*
    let mut edge_to_node: HashMap<EdgeIndex, NodeIndex> = HashMap::default();

    // add node for each edge
    for edge in graph.edge_indices() {
        // for a edge from node_a to node_b
        let (node_a, node_b) = graph.edge_endpoints(edge).unwrap();
        let node_weight_a = graph.node_weight(node_a).unwrap();
        let node_weight_b = graph.node_weight(node_b).unwrap();

        // add a node corresponds to the edge
        let node = g.add_node(node_map(node_weight_a, node_weight_b));
        edge_to_node.insert(edge, node);
    }

    // add edge for each node
    for node in graph.node_indices() {
        // for each pair of (incoming, outcoming) edges of a node
        //
        // G:
        //  edge_in        edge_out
        // --------> node --------->
        //
        //            |
        //            v
        // G*:
        //  node_in -----> node_out
        //
        for edge_in in graph.edges_directed(node, Direction::Incoming) {
            for edge_out in graph.edges_directed(node, Direction::Outgoing) {
                let node_in = *edge_to_node.get(&edge_in.id()).unwrap();
                let node_out = *edge_to_node.get(&edge_out.id()).unwrap();
                g.add_edge(node_in, node_out, ());
            }
        }
    }
    g
}

///
/// convert undirected graph into directed graph
/// by orienting an edge into higher indexed node.
///
fn to_directed<N: Clone, E: Clone>(graph: &UnGraph<N, E>) -> DiGraph<N, E> {
    let mut g = DiGraph::with_capacity(graph.node_count(), graph.edge_count());

    // add nodes
    for node in graph.node_indices() {
        let node_weight = graph.node_weight(node).unwrap().clone();
        g.add_node(node_weight);
    }

    // add edges
    for edge in graph.edge_indices() {
        let (v, w) = graph.edge_endpoints(edge).unwrap();
        let edge_weight = graph.edge_weight(edge).unwrap().clone();
        if v <= w {
            // v -> w
            g.add_edge(v, w, edge_weight);
        } else {
            // w -> v
            g.add_edge(w, v, edge_weight);
        }
    }

    g
}

///
/// create a edge list as a vector of (source, target) of the given directed graph.
///
fn to_edge_endpoint_list<N, E>(graph: &DiGraph<N, E>) -> Vec<(usize, usize)> {
    graph
        .edge_indices()
        .map(|e| {
            let (v, w) = graph.edge_endpoints(e).unwrap();
            (v.index(), w.index())
        })
        .collect()
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use petgraph::dot::Dot;

    #[test]
    fn cycle_space_into_adj_graph_01() {
        let g: UnGraph<(), ()> = UnGraph::from_edges(&[(0, 1), (2, 1), (2, 0)]);
        let g2 = to_directed(&g);
        println!("{:?}", Dot::with_config(&g2, &[]));
        assert_eq!(to_edge_endpoint_list(&g2), vec![(0, 1), (1, 2), (0, 2)]);
        let h = into_adj_graph(&g2, |_, _| ());
        println!("{:?}", Dot::with_config(&h, &[]));
        assert_eq!(to_edge_endpoint_list(&g2), vec![(0, 1)]);
    }

    #[test]
    fn cycle_space_to_directed() {
        let mut g: UnGraph<(), ()> = UnGraph::from_edges(&[(0, 1), (2, 1), (2, 0)]);
        println!("{:?}", Dot::with_config(&g, &[]));
        let g2 = to_directed(&g);
        println!("{:?}", Dot::with_config(&g2, &[]));
        assert_eq!(to_edge_endpoint_list(&g2), vec![(0, 1), (1, 2), (0, 2)]);
    }
}
