use fnv::FnvHashMap as HashMap;
use petgraph::algo::tarjan_scc;
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
#[derive(Clone, Debug)]
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
    ///
    /// Get composition mapping of two mappings (self and other)
    ///
    /// ```text
    /// f, g: EdgeMap
    /// A, B, C: DiGraph
    ///
    ///    f     g
    /// A --> B --> C
    /// ```
    ///
    pub fn compose(f: &EdgeMap, g: &EdgeMap) -> EdgeMap {
        let mut m = EdgeMap {
            edge_count_original: f.edge_count_original,
            edge_count: g.edge_count,
            from_original: HashMap::default(),
            to_original: HashMap::default(),
        };

        // from_original
        // (1) mappings in f
        for (&a, &b) in f.from_original.iter() {
            let c = b.and_then(|b| g.from_original(b));
            m.from_original.insert(a, c);
        }
        // (2) mappings not in f but in g
        for (&b, &c) in g.from_original.iter() {
            let a = f.to_original(b);
            m.from_original.insert(a, c);
        }

        // to_original
        // (1) mappings in g
        for (&c, &b) in g.to_original.iter() {
            let a = f.to_original(b);
            m.to_original.insert(c, a);
        }
        // (2) mappings not in g but in f
        for (&b, &a) in f.to_original.iter() {
            if let Some(c) = g.from_original(b) {
                m.to_original.insert(c, a);
            }
        }

        m
    }
}

/// Remove edges
///
///
pub fn purge_edges_with_mapping<N, E>(graph: &mut DiGraph<N, E>, edges: &[EdgeIndex]) -> EdgeMap {
    let mut edge_map = EdgeMap::new(graph);

    for &edge_remove_original in edges {
        edge_map.delete_original(graph, edge_remove_original);
    }

    edge_map
}

///
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
/// edges that bridges between two strongly connected components
///
pub fn bridge_edges<N, E>(graph: &DiGraph<N, E>) -> Vec<EdgeIndex> {
    // compute strongly connected components of graph by Tarjan's algorithm.
    // scc is a set of nodes in which any two pair of nodes v, w have a path v->w and w->v.
    let components = tarjan_scc(graph);

    // Convert
    //     components[id] = [node]
    // into
    //     component_ids[node] = id
    let mut component_ids: Vec<usize> = vec![0; graph.node_count()];
    for (id, component) in components.iter().enumerate() {
        for node in component {
            component_ids[node.index()] = id;
        }
    }

    // bridge edge is an edge that connects two nodes belonging to different component.
    //
    let mut ret = Vec::new();
    for edge in graph.edge_indices() {
        let (s, t) = graph.edge_endpoints(edge).unwrap();
        if component_ids[s.index()] != component_ids[t.index()] {
            ret.push(edge);
        }
    }

    ret
}

#[allow(dead_code)]
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
            let m = purge_edges_with_mapping(&mut g, &[ei(0)]);
            println!("{:?}", m);
            println!("{:?}", petgraph::dot::Dot::with_config(&g, &[]));
            assert_is_bimap(&m.from_original, &m.to_original);
            assert_eq!(
                m.from_original,
                HashMap::from_iter([(ei(0), None), (ei(8), Some(ei(0)))])
            );
            assert_eq!(m.to_original, HashMap::from_iter([(ei(0), ei(8))]));
            assert_eq!(g.edge_count(), 8);
            assert_eq!(g.node_count(), 8);
        }

        {
            let mut g = graph.clone();
            let m = purge_edges_with_mapping(&mut g, &[ei(3), ei(8), ei(5)]);
            println!("{:?}", m);
            println!("{:?}", petgraph::dot::Dot::with_config(&g, &[]));
            assert_is_bimap(&m.from_original, &m.to_original);
            assert_eq!(
                m.from_original,
                HashMap::from_iter([
                    (ei(3), None),
                    (ei(8), None),
                    (ei(5), None),
                    (ei(7), Some(ei(3))),
                    (ei(6), Some(ei(5))),
                ])
            );
            assert_eq!(
                m.to_original,
                HashMap::from_iter([(ei(5), ei(6)), (ei(3), ei(7))])
            );
            assert_eq!(g.edge_count(), 6);
            assert_eq!(g.node_count(), 8);
        }

        {
            let mut g = graph.clone();
            let m = purge_edges_with_mapping(
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
            println!("{:?}", m);
            println!("{:?}", petgraph::dot::Dot::with_config(&g, &[]));
            assert_is_bimap(&m.from_original, &m.to_original);
            assert_eq!(
                m.from_original,
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
            assert_eq!(m.to_original, HashMap::default());
            assert_eq!(g.edge_count(), 0);
            assert_eq!(g.node_count(), 8);
        }
    }
    #[test]
    fn edge_map_compose() {
        let mut graph: DiGraph<(), usize> = DiGraph::from_edges(&[
            // s, t, ei
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
        println!("a{:?}", petgraph::dot::Dot::with_config(&graph, &[]));
        let f = purge_edges_with_mapping(&mut graph, &[ei(3), ei(8), ei(5)]);
        println!("f{:?}", f);
        println!("b{:?}", petgraph::dot::Dot::with_config(&graph, &[]));
        let g = purge_edges_with_mapping(&mut graph, &[ei(1), ei(4)]);
        println!("c{:?}", petgraph::dot::Dot::with_config(&graph, &[]));
        println!("g{:?}", g);
        let h = EdgeMap::compose(&f, &g);
        println!("h{:?}", h);
        assert_is_bimap(&h.from_original, &h.to_original);
    }
}
