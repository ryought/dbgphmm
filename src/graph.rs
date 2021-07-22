#[derive(Debug, PartialEq, PartialOrd, Eq, Hash, Copy, Clone)]
pub struct Node(pub usize);

#[derive(Debug, PartialEq, PartialOrd, Eq, Hash, Copy, Clone)]
pub struct Edge(pub usize);

#[derive(Debug, Clone)]
pub struct IndexedDiGraph {
    n_nodes: usize,
    n_edges: usize,
    in_edges: Vec<Vec<Edge>>,
    out_edges: Vec<Vec<Edge>>,
    edges: Vec<(Node, Node)>,
}

impl IndexedDiGraph {
    fn from(edges: Vec<(Node, Node)>) -> IndexedDiGraph {
        // check the consistency
        let n_edges = edges.len();
        let n_nodes = edges.iter().map(|(v, w)| v.0.max(w.0)).max().unwrap() + 1;

        // create in_edges and out_edges
        let mut in_edges = vec![Vec::new(); n_nodes];
        let mut out_edges = vec![Vec::new(); n_nodes];
        for (i, (v, w)) in edges.iter().enumerate() {
            // e = (v, w): edge e from v to w
            let e = Edge(i);

            in_edges[w.0].push(e);
            out_edges[v.0].push(e);
        }

        IndexedDiGraph {
            n_edges,
            n_nodes,
            in_edges,
            out_edges,
            edges,
        }
    }
    fn n_nodes(&self) -> usize {
        self.n_nodes
    }
    fn n_edges(&self) -> usize {
        self.n_edges
    }
    fn in_edges(&self, v: &Node) -> &[Edge] {
        &self.in_edges[v.0]
    }
    fn out_edges(&self, v: &Node) -> &[Edge] {
        &self.out_edges[v.0]
    }
    fn node_pair(&self, e: &Edge) -> (Node, Node) {
        self.edges[e.0]
    }

    fn childs(&self, v: &Node) -> Vec<Node> {
        self.out_edges(v)
            .iter()
            .map(|e| {
                let (v2, w) = self.node_pair(e);
                assert_eq!(*v, v2);
                w
            })
            .collect()
    }
    fn parents(&self, v: &Node) -> Vec<Node> {
        self.in_edges(v)
            .iter()
            .map(|e| {
                let (w, v2) = self.node_pair(e);
                assert_eq!(*v, v2);
                w
            })
            .collect()
    }

    /// Floyd-Warshall-like dynamic programming to calculate
    /// F[k][v] = (minimum weight path from source to v, with k edges)
    /// with backtraces
    /// B[k][v] = w  <=>  (min weight path F[k][v] ends with w -> v edge)
    /// this function returns their products, that is
    /// R[k][v] = (F[k][v], B[k][v])
    fn min_weight_paths(&self, source: &Node, weights: &[f64]) -> Vec<Vec<(f64, Option<Node>)>> {
        let mut fs = Vec::new();
        // 1. initialize
        let mut f: Vec<(f64, Option<Node>)> = (0..self.n_nodes())
            .map(|i| {
                let v = Node(i);
                if v == *source {
                    (0.0, None)
                } else {
                    (f64::INFINITY, None)
                }
            })
            .collect();

        // 2. step
        for k in 1..self.n_nodes() {
            let f_new = (0..self.n_nodes())
                .map(|i| {
                    let v = Node(i);
                    if self.in_edges(&v).len() == 0 {
                        (f64::INFINITY, None)
                    } else {
                        self.in_edges(&v)
                            .iter()
                            .map(|e| {
                                let (w, _) = self.node_pair(e);
                                let weight = weights[e.0] + f[w.0].0;
                                if weight == f64::INFINITY {
                                    (weight, None)
                                } else {
                                    (weight, Some(w))
                                }
                            })
                            .min_by(|(a, _), (b, _)| a.partial_cmp(b).unwrap())
                            .unwrap()
                    }
                })
                .collect();
            fs.push(f);
            f = f_new;
        }

        // 3. finalize
        fs.push(f);
        fs
    }
    /*
    /// Find the min mean(average) weight cycle
    fn minimum_mean_weight_cycle(&self) -> Vec<Edge>;
    */
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn construction() {
        let v = vec![
            (Node(0), Node(1)),
            (Node(0), Node(2)),
            (Node(1), Node(3)),
            (Node(2), Node(3)),
            (Node(2), Node(4)),
        ];
        let g = IndexedDiGraph::from(v);

        assert_eq!(g.n_nodes(), 5);
        assert_eq!(g.n_edges(), 5);

        assert_eq!(g.in_edges(&Node(0)).len(), 0);
        assert_eq!(g.out_edges(&Node(0)).len(), 2);

        assert_eq!(g.parents(&Node(0)).len(), 0);
        assert_eq!(g.childs(&Node(0)).len(), 2);

        assert_eq!(g.parents(&Node(3)), vec![Node(1), Node(2)]);
        assert_eq!(g.childs(&Node(2)), vec![Node(3), Node(4)]);
        assert_eq!(g.childs(&Node(4)), vec![]);
    }

    #[test]
    fn min_weight_paths() {
        let v = vec![
            (Node(0), Node(1)),
            (Node(0), Node(2)),
            (Node(1), Node(3)),
            (Node(2), Node(3)),
            (Node(2), Node(4)),
        ];
        let g = IndexedDiGraph::from(v);
        let weights = vec![1.0, 1.0, 2.0, 1.0, 2.0];
        let paths = g.min_weight_paths(&Node(0), &weights);
        // println!("{:?}", paths);
        assert_eq!(paths[1][1], (weights[0], Some(Node(0)))); // 0->1
        assert_eq!(paths[1][2], (weights[1], Some(Node(0)))); // 0->2
        assert_eq!(paths[2][3], (weights[1] + weights[3], Some(Node(2)))); // 0->2->3
        assert_eq!(paths[2][4], (weights[1] + weights[4], Some(Node(2)))); // 0->2->4
    }
}
