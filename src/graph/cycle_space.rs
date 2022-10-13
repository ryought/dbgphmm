//!
//! cycle space
//!
//! enumerate all cycles in the undirected graph
//!
use super::cycle::{Cycle, SimpleCycle};
use fixedbitset::FixedBitSet;
use fnv::{FnvHashMap as HashMap, FnvHashSet as HashSet};
use itertools::Itertools;
use petgraph::graph::{DiGraph, EdgeIndex, Graph, NodeIndex, UnGraph};
use petgraph::stable_graph::IndexType;
use petgraph::visit::EdgeRef;
use petgraph::{Direction, EdgeType};

// unused imports
// use petgraph::unionfind::UnionFind;

#[derive(Clone, Debug)]
struct CycleSpace {
    ///
    /// a set of all cycle basis.
    ///
    basis: Vec<SimpleCycle>,
}

struct AdjGraphNodeIterator<N: PartialEq + Clone, F: Fn(&N, &N) -> N> {
    /// degree of adjacency (i.e. how many times the graph was upconverted)
    k: usize,
    /// node indices
    i: usize,
    /// adjacency graph
    graph: UnGraph<N, ()>,
    /// node mapping function
    node_map: F,
}

impl<N: PartialEq + Clone, F: Fn(&N, &N) -> N> AdjGraphNodeIterator<N, F> {
    /// constructor of `AdjGraphNodeIterator`
    /// start from k=1 and i=0.
    pub fn new(graph: UnGraph<N, ()>, node_map: F) -> Self {
        AdjGraphNodeIterator {
            k: 1,
            i: 0,
            graph,
            node_map,
        }
    }
    fn upconvert(&mut self) {
        self.graph = into_adj_graph_no_node_dups(&self.graph, |x, y| (self.node_map)(x, y));
        self.k += 1;
        self.i = 0;
    }
}

impl<N: PartialEq + Clone, F: Fn(&N, &N) -> N> Iterator for AdjGraphNodeIterator<N, F> {
    type Item = N;
    fn next(&mut self) -> Option<N> {
        // upconvert
        let n = self.graph.node_count();
        if n > 0 && self.i == n {
            self.upconvert()
        }

        // yield node[i]
        if self.graph.node_count() > 0 {
            let w = self.graph.node_weight(NodeIndex::new(self.i)).unwrap();
            self.i += 1;
            Some(w.clone())
        } else {
            None
        }
    }
}

///
/// convert a graph `G` into adjacent graph `G*`, with node duplication checks.
///
/// - an edge in `G` corresponds to a node in `G*`
/// - a pair of nodes (v, w) in `G*` will be connected by an edge v->w if corresponding edges in `G`
///   if the corresponding edge ends/starts from the same node.
///
fn into_adj_graph_no_node_dups<N, F>(graph: &UnGraph<N, ()>, node_map: F) -> UnGraph<N, ()>
where
    N: PartialEq,
    F: Fn(&N, &N) -> N,
{
    let mut g = UnGraph::new_undirected();
    // map from edge in G into node in G*
    let mut edge_to_node: HashMap<EdgeIndex, NodeIndex> = HashMap::default();

    // add node for each edge
    for edge in graph.edge_indices() {
        // for a edge from node_a to node_b
        let (node_a, node_b) = graph.edge_endpoints(edge).unwrap();
        let node_weight_a = graph.node_weight(node_a).unwrap();
        let node_weight_b = graph.node_weight(node_b).unwrap();

        // create a node corresponds to the edge if not exists
        let new_node_weight = node_map(node_weight_a, node_weight_b);

        let node = g.node_indices().find(|node| {
            let node_weight = g.node_weight(*node).unwrap();
            node_weight == &new_node_weight
        });

        match node {
            Some(node) => {
                // the node already exists
                edge_to_node.insert(edge, node);
            }
            None => {
                // the node is new and to be created
                let node = g.add_node(new_node_weight);
                edge_to_node.insert(edge, node);
            }
        }
    }

    // add edge for each node
    for node in graph.node_indices() {
        // for each pair of edges connected to the node
        for (e_a, e_b) in graph.edges(node).tuple_combinations() {
            let node_a = *edge_to_node.get(&e_a.id()).unwrap();
            let node_b = *edge_to_node.get(&e_b.id()).unwrap();
            if node_a != node_b && !g.contains_edge(node_a, node_b) {
                g.add_edge(node_a, node_b, ());
            }
        }
    }
    g
}

/// join a sorted list
fn test_node_map(xs: &[u8], ys: &[u8]) -> Vec<u8> {
    let mut ret = xs.to_vec();
    for y in ys {
        match ret.binary_search(y) {
            Ok(_) => {}                      // ignore the duplicated element
            Err(pos) => ret.insert(pos, *y), // add an element only in ys
        }
    }
    ret
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
fn to_edge_endpoint_list<N, E, Ty: EdgeType, Ix: IndexType>(
    graph: &Graph<N, E, Ty, Ix>,
) -> Vec<(usize, usize)> {
    graph
        .edge_indices()
        .map(|e| {
            let (v, w) = graph.edge_endpoints(e).unwrap();
            (v.index(), w.index())
        })
        .collect()
}

fn to_node_weight_list<N: Clone, E, Ty: EdgeType, Ix: IndexType>(
    graph: &Graph<N, E, Ty, Ix>,
) -> Vec<N> {
    graph
        .node_indices()
        .map(|n| graph.node_weight(n).unwrap().clone())
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
    fn cycle_space_into_adj_graph_no_node_dups_01() {
        // (1) test test_node_map
        let t1 = test_node_map(&[0, 1, 2], &[0, 2, 5]);
        println!("{:?}", t1);
        assert_eq!(t1, vec![0, 1, 2, 5]);

        // (2)
        let mut g: UnGraph<Vec<u8>, ()> = UnGraph::new_undirected();
        let v0 = g.add_node(vec![0]);
        let v1 = g.add_node(vec![1]);
        let v2 = g.add_node(vec![2]);
        g.add_edge(v0, v1, ());
        g.add_edge(v0, v2, ());
        g.add_edge(v1, v2, ());
        println!("{:?}", Dot::with_config(&g, &[]));
        assert_eq!(to_node_weight_list(&g), vec![vec![0], vec![1], vec![2]]);

        let h1 = into_adj_graph_no_node_dups(&g, |x, y| test_node_map(x, y));
        println!("{:?}", Dot::with_config(&h1, &[]));
        assert_eq!(
            to_node_weight_list(&h1),
            vec![vec![0, 1], vec![0, 2], vec![1, 2]]
        );

        let h2 = into_adj_graph_no_node_dups(&h1, |x, y| test_node_map(x, y));
        println!("{:?}", Dot::with_config(&h2, &[]));
        assert_eq!(to_node_weight_list(&h2), vec![vec![0, 1, 2]]);

        let h3 = into_adj_graph_no_node_dups(&h2, |x, y| test_node_map(x, y));
        println!("{:?}", Dot::with_config(&h3, &[]));
        assert_eq!(to_node_weight_list(&h3).len(), 0);

        // test iterator
        let agiter = AdjGraphNodeIterator::new(g, |x, y| test_node_map(x, y));
        let xs: Vec<Vec<u8>> = agiter.collect();
        println!("{:?}", xs);
        assert_eq!(
            xs,
            vec![
                vec![0],
                vec![1],
                vec![2],
                vec![0, 1],
                vec![0, 2],
                vec![1, 2],
                vec![0, 1, 2]
            ]
        );
    }

    #[test]
    fn cycle_space_into_adj_graph_no_node_dups_02() {
        let mut g: UnGraph<Vec<u8>, ()> = UnGraph::new_undirected();
        let v1 = g.add_node(vec![1]);
        let v2 = g.add_node(vec![2]);
        let v3 = g.add_node(vec![3]);
        let v4 = g.add_node(vec![4]);
        let v5 = g.add_node(vec![5]);
        g.add_edge(v1, v2, ());
        g.add_edge(v1, v3, ());
        g.add_edge(v2, v5, ());
        g.add_edge(v3, v4, ());
        g.add_edge(v3, v5, ());
        println!("{:?}", Dot::with_config(&g, &[]));
        assert_eq!(
            to_node_weight_list(&g),
            vec![vec![1], vec![2], vec![3], vec![4], vec![5]]
        );

        let h1 = into_adj_graph_no_node_dups(&g, |x, y| test_node_map(x, y));
        println!("{:?}", Dot::with_config(&h1, &[]));
        assert_eq!(
            to_node_weight_list(&h1),
            vec![vec![1, 2], vec![1, 3], vec![2, 5], vec![3, 4], vec![3, 5],]
        );

        let h2 = into_adj_graph_no_node_dups(&h1, |x, y| test_node_map(x, y));
        println!("{:?}", Dot::with_config(&h2, &[]));
        assert_eq!(
            to_node_weight_list(&h2),
            vec![
                vec![1, 2, 3],
                vec![1, 2, 5],
                vec![3, 4, 5],
                vec![1, 3, 5],
                vec![1, 3, 4],
                vec![2, 3, 5],
            ]
        );

        let h3 = into_adj_graph_no_node_dups(&h2, |x, y| test_node_map(x, y));
        println!("{:?}", Dot::with_config(&h3, &[]));
        assert_eq!(
            to_node_weight_list(&h3),
            vec![
                vec![1, 2, 3, 5],
                vec![1, 2, 3, 4],
                vec![1, 3, 4, 5],
                vec![2, 3, 4, 5],
            ]
        );

        let h4 = into_adj_graph_no_node_dups(&h3, |x, y| test_node_map(x, y));
        println!("{:?}", Dot::with_config(&h4, &[]));
        assert_eq!(to_node_weight_list(&h4), vec![vec![1, 2, 3, 4, 5]]);

        // test iterator
        let agiter = AdjGraphNodeIterator::new(g, |x, y| test_node_map(x, y));
        let xs: Vec<Vec<u8>> = agiter.collect();
        assert_eq!(
            xs,
            vec![
                vec![1],
                vec![2],
                vec![3],
                vec![4],
                vec![5],
                vec![1, 2],
                vec![1, 3],
                vec![2, 5],
                vec![3, 4],
                vec![3, 5],
                vec![1, 2, 3],
                vec![1, 2, 5],
                vec![3, 4, 5],
                vec![1, 3, 5],
                vec![1, 3, 4],
                vec![2, 3, 5],
                vec![1, 2, 3, 5],
                vec![1, 2, 3, 4],
                vec![1, 3, 4, 5],
                vec![2, 3, 4, 5],
                vec![1, 2, 3, 4, 5]
            ]
        );
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
