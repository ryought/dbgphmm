use fnv::FnvHashMap as HashMap;
use petgraph::graph::{DiGraph, EdgeIndex, Graph, NodeIndex, UnGraph};
use petgraph::visit::EdgeRef;
use petgraph::Direction;
use std::iter::FromIterator;

pub fn to_node_list<N: Clone, E>(graph: &DiGraph<N, E>) -> Vec<(usize, N)> {
    graph
        .node_indices()
        .map(|node| {
            let weight = graph.node_weight(node).unwrap().clone();
            (node.index(), weight)
        })
        .collect()
}

pub fn to_edge_list<N, E: Clone>(graph: &DiGraph<N, E>) -> Vec<(usize, usize, usize, E)> {
    graph
        .edge_indices()
        .map(|edge| {
            let weight = graph.edge_weight(edge).unwrap().clone();
            let (source, target) = graph.edge_endpoints(edge).unwrap();
            (edge.index(), source.index(), target.index(), weight)
        })
        .collect()
}

///
/// count the number of (in_degree, out_degree) of nodes.
///
pub fn degree_stats<N, E>(graph: &DiGraph<N, E>) -> HashMap<(usize, usize), usize> {
    let mut h = HashMap::default();
    for node in graph.node_indices() {
        let in_degree = graph.edges_directed(node, Direction::Incoming).count();
        let out_degree = graph.edges_directed(node, Direction::Outgoing).count();
        *h.entry((in_degree, out_degree)).or_insert(0) += 1;
    }
    h
}

///
/// Store correspondence of edges of the graph before and after edge deletions
///
pub struct EdgeMap {
    ///
    ///
    ///
    edge_count: usize,
    ///
    ///
    ///
    edge_count_original: usize,
    ///
    /// mapping of edge in original graph -> edge in current graph
    ///
    /// from_original[edge in original] = None  <=>  edge in original is deleted
    ///
    from_original: HashMap<EdgeIndex, Option<EdgeIndex>>,
    ///
    /// mapping of edge in current graph (after removing) -> edge in original graph (before
    /// removing)
    ///
    /// to_original[edge in current] = edge in original
    ///
    to_original: HashMap<EdgeIndex, EdgeIndex>,
}

impl EdgeMap {
    ///
    ///
    pub fn new<N, E>(graph: &DiGraph<N, E>) -> EdgeMap {
        EdgeMap {
            edge_count: graph.edge_count(),
            edge_count_original: graph.edge_count(),
            from_original: HashMap::default(),
            to_original: HashMap::default(),
        }
    }
    ///
    ///
    ///
    pub fn delete<N, E>(&mut self, graph: &mut DiGraph<N, E>, edge: EdgeIndex) {
        assert!(edge.index() < self.edge_count, "edge is out of index");

        // (A) update EdgeMap
        let edge_last = EdgeIndex::new(graph.edge_count() - 1);
        if edge == edge_last {
            // the edge to be removed is the last edge in the graph
            // index of other edges left unmodified

            // (1) update to_original
            let edge_original = self.to_original(edge);
            self.to_original.remove(&edge);

            // (2) update from_original
            self.from_original.insert(edge_original, None);
        } else {
            // the edge to be removed is not the last edge in the graph
            // so the last edge will be
            let edge_last_original = self.to_original(edge_last);
            let edge_original = self.to_original(edge);

            // (1) update to_original
            self.to_original.insert(edge, edge_last_original);
            self.to_original.remove(&edge_last);

            // (2) update from_original
            self.from_original.insert(edge_original, None);
            self.from_original.insert(edge_last_original, Some(edge));
        }

        // (B) update DiGraph
        graph
            .remove_edge(edge)
            .expect("edge to be removed does not exist");

        // (C) update edge_count
        self.edge_count -= 1;
    }
    ///
    ///
    ///
    pub fn delete_original<N, E>(&mut self, graph: &mut DiGraph<N, E>, edge_original: EdgeIndex) {
        assert!(
            edge_original.index() < self.edge_count_original,
            "edge is out of index, {} {}",
            edge_original.index(),
            self.edge_count_original,
        );

        if let Some(edge) = self.from_original(edge_original) {
            self.delete(graph, edge);
        }
    }
    ///
    /// map edge in current graph into edge in original graph
    ///
    pub fn to_original(&self, edge: EdgeIndex) -> EdgeIndex {
        assert!(edge.index() < self.edge_count, "edge is out of index");

        self.to_original.get(&edge).copied().unwrap_or(edge)
    }
    ///
    /// map edge in original graph into edge in current graph
    ///
    pub fn from_original(&self, edge: EdgeIndex) -> Option<EdgeIndex> {
        assert!(
            edge.index() < self.edge_count_original,
            "edge is out of index"
        );

        match self.from_original.get(&edge) {
            None => {
                // edge is unmodified after the process
                Some(edge)
            }
            Some(None) => {
                // edge is deleted
                None
            }
            Some(Some(edge)) => {
                // edge index is changed
                Some(*edge)
            }
        }
    }
}

/// Remove edges
///
///
pub fn purge_edges_with_mapping<N, E>(
    graph: &mut DiGraph<N, E>,
    edges: &[EdgeIndex],
) -> (
    HashMap<EdgeIndex, Option<EdgeIndex>>,
    HashMap<EdgeIndex, EdgeIndex>,
) {
    let mut edge_map = EdgeMap::new(graph);

    for &edge_remove_original in edges {
        edge_map.delete_original(graph, edge_remove_original);
    }

    (edge_map.from_original, edge_map.to_original)
}

/// Delete all isolated nodes (no in/out edges)
///
pub fn delete_isolated_nodes<N, E>(graph: &mut DiGraph<N, E>) {
    graph.retain_nodes(|graph, node| {
        // keep nodes with in/out edges, and remove otherwise
        let in_degree = graph.edges_directed(node, Direction::Incoming).count();
        let out_degree = graph.edges_directed(node, Direction::Outgoing).count();
        in_degree > 0 || out_degree > 0
    })
}

///
///
///
pub fn delete_unreachable_edges<N, E>(graph: &mut DiGraph<N, E>) {
    // delete nodes with no in/out edges
    unimplemented!();
}

fn assert_is_bimap(
    from_original: &HashMap<EdgeIndex, Option<EdgeIndex>>,
    to_original: &HashMap<EdgeIndex, EdgeIndex>,
) {
    for (x, y) in to_original.iter() {
        assert_eq!(from_original.get(y).unwrap().unwrap(), *x);
    }
}

///
/// Split a node into two nodes with preserving edge ids
///
/// ```text
/// v1 -->      --> w1
///        node
/// v2 -->      --> w2
///
/// into
///
/// v1 -->                      --> w1
///        node_in --> node_out
/// v2 -->                      --> w2
/// ```
///
pub fn split_node<N: Clone, E: Clone>(graph: &mut DiGraph<N, E>, node: NodeIndex, edge_weight: E) {
    let node_weight = graph.node_weight(node).unwrap().clone();

    // add two new nodes and an edge between them
    let node_in = graph.add_node(node_weight.clone());
    let node_out = graph.add_node(node_weight.clone());
    graph.add_edge(node_in, node_out, edge_weight);

    // connect in-edges of node
    //   v --> node
    // into
    //   v --> node_in -> node_out
    let edges_in: Vec<_> = graph
        .edges_directed(node, Direction::Incoming)
        .map(|e| e.id())
        .collect();
    for edge_in in edges_in {
        let (source, _) = graph.edge_endpoints(edge_in).unwrap();
        let weight = graph.edge_weight(edge_in).unwrap().clone();
        // change edge while preserving id
        graph.add_edge(source, node_in, weight);
        graph.remove_edge(edge_in);
    }

    // connect out-edges of node
    //   node --> w
    // into
    //   node_in -> node_out -> w
    let edges_out: Vec<_> = graph
        .edges_directed(node, Direction::Outgoing)
        .map(|e| e.id())
        .collect();
    for edge_out in edges_out {
        let (_, target) = graph.edge_endpoints(edge_out).unwrap();
        let weight = graph.edge_weight(edge_out).unwrap().clone();
        // change edge while preserving id
        graph.add_edge(node_out, target, weight);
        graph.remove_edge(edge_out);
    }

    // remove the splitted node
    graph.remove_node(node);
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ei;

    #[test]
    fn purge_edges_with_mapping_test() {
        let graph: DiGraph<(), usize> = DiGraph::from_edges(&[
            (0, 1, 0),
            (1, 2, 1),
            (2, 3, 2),
            (4, 3, 3),
            (2, 5, 4),
            (3, 7, 5),
            (5, 7, 6),
            (5, 6, 7),
            (6, 5, 8),
        ]);
        assert_eq!(graph.edge_count(), 9);
        assert_eq!(graph.node_count(), 8);

        {
            let mut g = graph.clone();
            let (f, t) = purge_edges_with_mapping(&mut g, &[ei(0)]);
            println!("from_original={:?} to_original={:?}", f, t);
            println!("{:?}", petgraph::dot::Dot::with_config(&g, &[]));
            assert_is_bimap(&f, &t);
            assert_eq!(f, HashMap::from_iter([(ei(0), None), (ei(8), Some(ei(0)))]));
            assert_eq!(t, HashMap::from_iter([(ei(0), ei(8))]));
            assert_eq!(g.edge_count(), 8);
            assert_eq!(g.node_count(), 8);
        }

        {
            let mut g = graph.clone();
            let (f, t) = purge_edges_with_mapping(&mut g, &[ei(3), ei(8), ei(5)]);
            println!("from_original={:?} to_original={:?}", f, t);
            println!("{:?}", petgraph::dot::Dot::with_config(&g, &[]));
            assert_is_bimap(&f, &t);
            assert_eq!(
                f,
                HashMap::from_iter([
                    (ei(3), None),
                    (ei(8), None),
                    (ei(5), None),
                    (ei(7), Some(ei(3))),
                    (ei(6), Some(ei(5))),
                ])
            );
            assert_eq!(t, HashMap::from_iter([(ei(5), ei(6)), (ei(3), ei(7))]));
            assert_eq!(g.edge_count(), 6);
            assert_eq!(g.node_count(), 8);
        }

        {
            let mut g = graph.clone();
            let (f, t) = purge_edges_with_mapping(
                &mut g,
                &[
                    ei(3),
                    ei(8),
                    ei(5),
                    ei(7),
                    ei(0),
                    ei(4),
                    ei(6),
                    ei(1),
                    ei(2),
                ],
            );
            println!("from_original={:?} to_original={:?}", f, t);
            println!("{:?}", petgraph::dot::Dot::with_config(&g, &[]));
            assert_is_bimap(&f, &t);
            assert_eq!(
                f,
                HashMap::from_iter([
                    (ei(0), None),
                    (ei(1), None),
                    (ei(2), None),
                    (ei(3), None),
                    (ei(4), None),
                    (ei(5), None),
                    (ei(6), None),
                    (ei(7), None),
                    (ei(8), None),
                ])
            );
            assert_eq!(t, HashMap::default());
            assert_eq!(g.edge_count(), 0);
            assert_eq!(g.node_count(), 8);
        }
    }
}
