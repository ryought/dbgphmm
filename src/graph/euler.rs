//!
//! Calculating the number of euler circuits
//!
use crate::utils::log_factorial;
use ndarray::prelude::*;
use ndarray_linalg::solve::Determinant;
use petgraph::algo::connected_components;
use petgraph::algo::tarjan_scc;
use petgraph::graph::{DiGraph, EdgeIndex, Graph, NodeIndex, UnGraph};
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
///
///
pub fn euler_circuit_count(graph: &DiGraph<(), usize>) -> f64 {
    let mut graph = graph.clone();
    graph.retain_edges(|g, e| g[e] > 0);
    graph.retain_nodes(|g, v| g.edges_directed(v, Direction::Outgoing).count() > 0);

    println!("cc={}", connected_components(&graph));
    println!("cc={:?}", tarjan_scc(&graph));

    let mut ret = 0.0;
    let n = graph.node_count();
    if n == 0 {
        return f64::NEG_INFINITY;
    }

    // for each strongly connected component..
    for component in tarjan_scc(&graph) {
        let h = subgraph(&graph, &component);
        let count = euler_circuit_count_in_connected(&h);
        println!("compont={:?} count={}", component, count.exp());
        ret += count;
    }

    ret
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

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
}
