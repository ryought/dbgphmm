//!
//! Calculating the number of euler circuits
//!
use crate::utils::log_factorial;
use fnv::FnvHashMap as HashMap;
use ndarray::prelude::*;
use ndarray_linalg::solve::Determinant;
use petgraph::algo::tarjan_scc;
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
use petgraph::visit::EdgeRef;
use petgraph::Direction;

///
///
///
pub fn subgraph<N: Clone, E: Clone>(graph: &DiGraph<N, E>, nodes: &[NodeIndex]) -> DiGraph<N, E> {
    let mut graph = graph.clone();
    graph.retain_nodes(|_, node| nodes.contains(&node));
    graph
}

pub fn euler_circuit_count_in_connected(graph: &DiGraph<(), usize>) -> f64 {
    // println!("[euler] {:?}", petgraph::dot::Dot::with_config(&graph, &[]));
    let n = graph.node_count();
    if n == 0 {
        return 0.0;
    }

    //
    // PartA: create laplacian matrix L
    //
    let mut laplacian: Array2<f64> = Array::zeros((n, n));
    // (1) diag (degree) matrix
    // L[i,i] += (total copy numbers of node i)
    for i in graph.node_indices() {
        let c: usize = graph
            .edges_directed(i, Direction::Outgoing)
            .map(|e| e.weight())
            .sum();
        // println!("i={} c={}", i.index(), c);
        laplacian[[i.index(), i.index()]] = c as f64;
    }
    // (2) subtract adjacency matrix a_ij
    // L[i,j] -= (total copy numbers of edges i->j)
    for i in graph.node_indices() {
        for j in graph.node_indices() {
            let c: usize = graph.edges_connecting(i, j).map(|e| e.weight()).sum();
            // println!("i={} j={} c={}", i.index(), j.index(), c);
            laplacian[[i.index(), j.index()]] -= c as f64;
        }
    }
    // add +1
    // starting point is arbitrary
    laplacian[[0, 0]] += 1.0;
    // println!("L={}", laplacian);
    let (sign, ln) = laplacian.sln_det().unwrap();
    // println!("{} {}", sign, ln);
    // println!("detL={}", sign * ln.exp());

    //
    // PartB:
    //
    let mut count = if ln == f64::NEG_INFINITY {
        0.0
    } else {
        sign * ln
    };
    for i in graph.node_indices() {
        let c: usize = graph
            .edges_directed(i, Direction::Outgoing)
            .map(|e| e.weight())
            .sum();
        if c > 0 {
            count += log_factorial(c - 1);

            for e in graph.edges_directed(i, Direction::Outgoing) {
                let c = e.weight();
                count -= log_factorial(*c);
            }
        }
    }

    count
}

///
/// Count the number of Eulerian circuits in the graph (whose edge has its multiplicity) (in log
/// space)
///
/// if `allow_multiple_component=true`, the final count will be the product of count of each components.
/// otherwise, the resulting count is log(0) if the graph is not strongly connected. (because there
/// is no Eulerian circuit)
///
pub fn euler_circuit_count(graph: &DiGraph<(), usize>, allow_multiple_component: bool) -> f64 {
    let mut graph = graph.clone();
    // remove zero edges
    graph.retain_edges(|g, e| g[e] > 0);
    // remove isolated nodes (= no outgoing edge)
    graph.retain_nodes(|g, v| g.edges_directed(v, Direction::Outgoing).count() > 0);

    // println!("cc={}", connected_components(&graph));
    // println!("cc={:?}", tarjan_scc(&graph));

    let mut ret = 0.0;
    let n = graph.node_count();
    if n == 0 {
        return f64::NEG_INFINITY;
    }

    // for each strongly connected component..
    if allow_multiple_component {
        for component in tarjan_scc(&graph) {
            let h = subgraph(&graph, &component);
            let count = euler_circuit_count_in_connected(&h);
            // println!("compont={:?} count={}", component, count.exp());
            ret += count;
        }
    } else {
        if tarjan_scc(&graph).len() > 1 {
            // if there are separate component, it has no euler circuit
            return f64::NEG_INFINITY;
        } else {
            ret = euler_circuit_count_in_connected(&graph);
        }
    }

    ret
}

///
/// join two circuit at an intersection point
///
pub fn join_circuit(
    graph: &DiGraph<(), usize>,
    mut a: Vec<EdgeIndex>,
    mut b: Vec<EdgeIndex>,
) -> Vec<EdgeIndex> {
    // find intersection point (node that is passed by both cycle)
    let hash_a: HashMap<NodeIndex, usize> = a
        .iter()
        .enumerate()
        .map(|(i, &ea)| {
            let (source, _) = graph.edge_endpoints(ea).unwrap();
            (source, i)
        })
        .collect();
    let (ia, ib) = b
        .iter()
        .enumerate()
        .find_map(|(ib, &eb)| {
            let (source, _) = graph.edge_endpoints(eb).unwrap();
            hash_a.get(&source).map(|&ia| (ia, ib))
        })
        .expect("two cycle has no intersection");

    // a
    //
    // a        a2
    // 0 1 2    3 4
    // ->->-> v ->->
    // a  = [0, ia)
    // a2 = [ia, na)
    let mut a2 = a.split_off(ia);
    // b
    //
    // b        b2
    // 0 1 2    3 4
    // ->->-> v ->->
    // b  = [0, ib)
    // b2 = [ib, nb)
    let mut b2 = b.split_off(ib);

    // join
    // a     b2   b   a
    // ----> ---> --> ------->
    a.append(&mut b2);
    a.append(&mut b);
    a.append(&mut a2);

    a
}

///
/// get an Euler circuit as edge vector
///
pub fn euler_circuit(graph: &DiGraph<(), usize>, from: NodeIndex) -> Vec<EdgeIndex> {
    // remain[edge] = (how many multiplicity remains in the edge?)
    let mut remain: Vec<usize> = graph.edge_indices().map(|e| graph[e]).collect();
    let mut total_remain: usize = remain.iter().sum();
    let mut cycles: Vec<Vec<EdgeIndex>> = Vec::new();
    let node_remain = |node: NodeIndex, remain: &[usize]| {
        graph
            .edges_directed(node, Direction::Outgoing)
            .map(|edge| remain[edge.id().index()])
            .sum::<usize>()
    };

    // empty cycle is fine
    if total_remain == 0 {
        return vec![];
    }

    assert!(
        node_remain(from, &remain) > 0,
        "from node do not have outgoing edges"
    );

    let mut start = from;
    loop {
        // traverse to find a cycle
        let mut cycle: Vec<EdgeIndex> = Vec::new();
        let mut node = start;
        loop {
            match graph
                .edges_directed(node, Direction::Outgoing)
                .find(|edge| remain[edge.id().index()] > 0)
            {
                Some(edge) => {
                    cycle.push(edge.id());
                    remain[edge.id().index()] -= 1;
                    total_remain -= 1;
                    node = edge.target();
                }
                None => break,
            }
        }
        assert_eq!(
            start, node,
            "end node of cycle is different from start node"
        );
        cycles.push(cycle);

        // there is no more cycles
        if total_remain == 0 {
            break;
        }

        // find a new start node
        start = graph
            .node_indices()
            .find(|&node| node_remain(node, &remain) > 0)
            .unwrap();
    }

    let cycle = cycles
        .into_iter()
        .reduce(|cycle_a, cycle_b| join_circuit(graph, cycle_a, cycle_b))
        .unwrap();

    // check that cycle is passes all edges
    let total: usize = graph.edge_indices().map(|e| graph[e]).sum();
    assert_eq!(total, cycle.len());

    cycle
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ei;

    #[test]
    fn subgraph_test() {
        let mut g = DiGraph::new();
        let a = g.add_node(());
        let b = g.add_node(());
        let c = g.add_node(());
        g.add_edge(a, b, 1);
        g.add_edge(a, c, 2);
        g.add_edge(c, b, 3);
        g.add_edge(b, c, 4);
        g.add_edge(b, c, 5);
        let h = subgraph(&g, &[b, c]);
        println!("{:?}", petgraph::dot::Dot::with_config(&h, &[]));
    }

    #[test]
    fn n_euler() {
        let assert_euler_approx_eq = |graph: &DiGraph<(), usize>, n0: usize, n1: usize| {
            assert!((euler_circuit_count(&graph, false).exp() - n0 as f64).abs() < 0.001);
            assert!((euler_circuit_count(&graph, true).exp() - n1 as f64).abs() < 0.001);
        };

        // self loop
        let g: DiGraph<(), usize> = DiGraph::from_edges(&[(0, 0, 1)]);
        assert_euler_approx_eq(&g, 1, 1);

        // single edge has no euler circuit
        let g: DiGraph<(), usize> = DiGraph::from_edges(&[(0, 1, 2)]);
        assert_euler_approx_eq(&g, 0, 0);

        // loop
        let g: DiGraph<(), usize> = DiGraph::from_edges(&[(0, 1, 1), (1, 0, 1)]);
        assert_euler_approx_eq(&g, 1, 1);

        // three euler circuit
        // - A,X,X,B
        // - A,X,B,X
        // - A,B,X,X
        let g: DiGraph<(), usize> = DiGraph::from_edges(&[
            (0, 0, 1), // A
            (0, 0, 1), // B
            (0, 0, 2), // X
        ]);
        assert_euler_approx_eq(&g, 3, 3);

        // two bubbles with two euler circuit
        // - A1 B C1 D A2 B C2 D
        // - A1 B C2 D A2 B C1 D
        let g: DiGraph<(), usize> = DiGraph::from_edges(&[
            (0, 1, 1), // A1
            (0, 1, 1), // A2
            (1, 2, 2), // B
            (2, 3, 1), // C1
            (2, 3, 1), // C2
            (3, 0, 2), // D
        ]);
        assert_euler_approx_eq(&g, 2, 2);

        // separate components
        // if allowing separate components, it has euler circuit A and B
        // otherwise it has no circuit
        let g: DiGraph<(), usize> = DiGraph::from_edges(&[
            (0, 0, 1), // A
            (1, 1, 1), // B
        ]);
        assert_euler_approx_eq(&g, 0, 1);

        println!(
            "{} {}",
            euler_circuit_count(&g, false).exp(),
            euler_circuit_count(&g, true).exp(),
        );
    }

    #[test]
    fn euler_circuit_test() {
        // two bubbles with two euler circuit
        // - A1 B C1 D A2 B C2 D
        // - A1 B C2 D A2 B C1 D
        let g: DiGraph<(), usize> = DiGraph::from_edges(&[
            (0, 1, 1), // 0 A1
            (0, 1, 1), // 1 A2
            (1, 2, 2), // 2 B
            (2, 3, 1), // 3 C1
            (2, 3, 1), // 4 C2
            (3, 4, 1), // 5 L
            (4, 5, 1), // 6 L
            (5, 3, 1), // 7 L
            (3, 0, 2), // 8 D
        ]);

        let c = join_circuit(
            &g,
            vec![ei(0), ei(2), ei(4), ei(8)],
            vec![ei(6), ei(7), ei(5)],
        );
        assert_eq!(c, vec![ei(0), ei(2), ei(4), ei(5), ei(6), ei(7), ei(8)]);

        let c = join_circuit(
            &g,
            vec![ei(8), ei(0), ei(2), ei(4)],
            vec![ei(5), ei(6), ei(7)],
        );
        assert_eq!(c, vec![ei(5), ei(6), ei(7), ei(8), ei(0), ei(2), ei(4)]);

        let c = euler_circuit(&g, NodeIndex::new(0));
        assert_eq!(
            c,
            vec![
                ei(1),
                ei(2),
                ei(4),
                ei(8),
                ei(0),
                ei(2),
                ei(3),
                ei(5),
                ei(6),
                ei(7),
                ei(8)
            ]
        );
    }
}
