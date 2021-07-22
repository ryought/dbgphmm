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

    /*
    /// F[k][v] = ()
    fn shortest_path_weights(&self, source: &Node, weights: Vec<f64>) -> Vec<Vec<f64>>;
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
        // println!("{:#?}", g);

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
}
