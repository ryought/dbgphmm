use forceatlas2;
use log::info;
use serde::Serialize;

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

#[derive(Debug, Copy, Clone, Serialize)]
pub struct Pos(pub f32, pub f32);

impl IndexedDiGraph {
    pub fn from(edges: Vec<(Node, Node)>) -> IndexedDiGraph {
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
    pub fn n_nodes(&self) -> usize {
        self.n_nodes
    }
    pub fn n_edges(&self) -> usize {
        self.n_edges
    }
    pub fn in_edges(&self, v: &Node) -> &[Edge] {
        &self.in_edges[v.0]
    }
    pub fn out_edges(&self, v: &Node) -> &[Edge] {
        &self.out_edges[v.0]
    }
    pub fn node_pair(&self, e: &Edge) -> (Node, Node) {
        self.edges[e.0]
    }

    pub fn childs(&self, v: &Node) -> Vec<Node> {
        self.out_edges(v)
            .iter()
            .map(|e| {
                let (v2, w) = self.node_pair(e);
                assert_eq!(*v, v2);
                w
            })
            .collect()
    }
    pub fn parents(&self, v: &Node) -> Vec<Node> {
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
    /// for all 0<=k<=n and v \in V
    /// this function returns their products, that is
    /// R[k][v] = (F[k][v], B[k][v])
    pub fn min_weight_paths(
        &self,
        source: &Node,
        weights: &[f64],
    ) -> Vec<Vec<(f64, Option<(Node, Edge)>)>> {
        let mut fs = Vec::new();
        // 1. initialize
        let mut f: Vec<(f64, Option<(Node, Edge)>)> = (0..self.n_nodes())
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
        for k in 1..=self.n_nodes() {
            let f_new = (0..self.n_nodes())
                .map(|i| {
                    let v = Node(i);
                    if self.in_edges(&v).len() == 0 {
                        (f64::INFINITY, None)
                    } else {
                        self.in_edges(&v)
                            .iter()
                            .map(|&e| {
                                let (w, _) = self.node_pair(&e);
                                let weight = weights[e.0] + f[w.0].0;
                                if weight == f64::INFINITY {
                                    (weight, None)
                                } else {
                                    (weight, Some((w, e)))
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

    /// compute minimizer pair `(k, v)` that satisfies
    /// `(v*, k*) = argmin_v argmax_k (F[n][v] - F[k][v]) / (n-k)`
    /// it fails when there are no cycles in the graph.
    fn find_minimizer_pair(
        &self,
        source: &Node,
        weights: &[f64],
        paths: &[Vec<(f64, Option<(Node, Edge)>)>],
    ) -> Option<(Node, usize, f64)> {
        let n = self.n_nodes();
        let maxs: Vec<Option<(usize, f64)>> = (0..n)
            .map(|i| {
                (0..n)
                    .filter_map(|k| {
                        let fnv = paths[n][i].0;
                        let fkv = paths[k][i].0;
                        if fnv != f64::INFINITY && fkv != f64::INFINITY {
                            Some((k, (fnv - fkv) / (n - k) as f64))
                        } else {
                            None
                        }
                    })
                    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            })
            .collect();

        maxs.iter()
            .enumerate()
            .filter_map(|(i, max)| {
                if max.is_some() {
                    let (k, weight) = max.unwrap();
                    Some((Node(i), k, weight))
                } else {
                    None
                }
            })
            .min_by(|(_, _, a), (_, _, b)| a.partial_cmp(b).unwrap())
    }

    fn traceback_paths(
        &self,
        start: &Node,
        weights: &[f64],
        paths: &[Vec<(f64, Option<(Node, Edge)>)>],
    ) -> Vec<Edge> {
        // now: pointer to node
        let mut now = *start;
        // visited: (visited[node] = [i0, i1, ...]) means ()
        let mut visited: Vec<Option<usize>> = vec![None; self.n_nodes()];
        let mut edges: Vec<Edge> = Vec::new();

        for i in (0..paths.len()).rev() {
            // now in i-th node in cycle

            // 1. compare with previous occurrences j
            if let Some(j) = visited[now.0] {
                // ((paths.len() - 1) - j)-th element in nodes
                // corresponds to paths[j]
                let cycle = edges[(paths.len() - 1) - j..]
                    .iter()
                    .rev()
                    .map(|&v| v)
                    .collect();
                return cycle;
            }

            // 2. mark as visited
            visited[now.0] = Some(i);

            // 3. traceback
            if i != 0 {
                let (_, parent) = paths[i][now.0];
                now = parent.unwrap().0;
                edges.push(parent.unwrap().1);
            }
        }

        panic!("traceback failed");
    }

    pub fn cycle_weight(&self, cycle: &[Edge], weights: &[f64]) -> f64 {
        cycle.iter().map(|e| weights[e.0]).sum()
    }

    pub fn cycle_as_node_list(&self, cycle: &[Edge]) -> Vec<Node> {
        cycle.iter().map(|e| self.node_pair(e).0).collect()
    }

    /// Find the min mean(average) weight cycle
    pub fn minimum_mean_weight_cycle(&self, source: &Node, weights: &[f64]) -> Option<Vec<Edge>> {
        // 1. compute shortest-paths
        let paths = self.min_weight_paths(source, weights);

        // 2. find minimizer
        let min = self.find_minimizer_pair(source, weights, &paths);

        // 3. trackback to find the min-mean-weight cycle
        match min {
            Some((v, k, expected_mean_weight)) => {
                let cycle = self.traceback_paths(&v, weights, &paths);
                let mean_weight = self.cycle_weight(&cycle, weights) / cycle.len() as f64;

                // check if this cycle has the desired mean_weight
                // assert_abs_diff_eq!(mean_weight, expected_mean_weight);

                Some(cycle)
            }
            None => {
                // no cycles
                None
            }
        }
    }

    /// calculate 2d layout
    pub fn layout_by_force_atlas2(&self) -> Vec<Pos> {
        const ITERATIONS: u32 = 5000;
        let edges: Vec<(usize, usize)> = self.edges.iter().map(|(v, w)| (v.0, w.0)).collect();
        info!("n_edges={}", edges.len());
        let mut layout = forceatlas2::Layout::<f32>::from_graph(
            edges,
            forceatlas2::Nodes::Degree(self.n_nodes),
            None,
            forceatlas2::Settings {
                chunk_size: None,
                dimensions: 2,
                dissuade_hubs: false,
                ka: 0.01,
                kg: 0.001,
                kr: 0.002,
                lin_log: true,
                speed: 1.0,
                prevent_overlapping: None,
                strong_gravity: false,
            },
        );
        info!("determining layout by forceatlas2...");
        for _ in 0..ITERATIONS {
            layout.iteration();
        }
        info!("n_nodes={}", layout.points.points.len());
        layout
            .points
            .iter()
            .map(|pos| Pos(pos[0], pos[1]))
            .collect()
    }
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

    /// Example graph 1
    /// 0 ---> 1 --+
    /// |          |
    /// +--> 2 ----+--> 3
    ///      |
    ///      +--> 4
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
        assert_eq!(paths[1][1], (weights[0], Some((Node(0), Edge(0))))); // 0->1
        assert_eq!(paths[1][2], (weights[1], Some((Node(0), Edge(1))))); // 0->2
        assert_eq!(
            paths[2][3],
            (weights[1] + weights[3], Some((Node(2), Edge(3))))
        ); // 0->2->3
        assert_eq!(
            paths[2][4],
            (weights[1] + weights[4], Some((Node(2), Edge(4))))
        ); // 0->2->4
    }

    /// counterexample in 'A note on finding minimum mean cycle' (Chaturvedi, 2017)
    #[test]
    fn min_mean_weight_cycle_01() {
        let v = vec![
            (Node(0), Node(1)),
            (Node(1), Node(2)),
            (Node(2), Node(0)),
            (Node(1), Node(3)),
            (Node(3), Node(4)),
            (Node(4), Node(5)),
            (Node(5), Node(6)),
            (Node(6), Node(1)),
            (Node(3), Node(7)),
            (Node(7), Node(4)),
        ];
        let g = IndexedDiGraph::from(v);
        let weights = vec![1.0, 3.0, -1.0, 2.0, 1.0, -1.0, 2.0, 1.0, 1.0, 2.0];
        let cycle = g.minimum_mean_weight_cycle(&Node(0), &weights).unwrap();
        assert_eq!(cycle, vec![Edge(3), Edge(4), Edge(5), Edge(6), Edge(7)]);
    }

    /// more simple graph with two cycles
    ///         3 ----> 4
    ///         ^       |
    ///         |       V
    /// 0 ----> 1 <---- 5
    /// ^       |
    /// +-- 2 <-+
    #[test]
    fn min_mean_weight_cycle_02() {
        let v = vec![
            (Node(0), Node(1)),
            (Node(1), Node(2)),
            (Node(2), Node(0)),
            (Node(1), Node(3)),
            (Node(3), Node(4)),
            (Node(4), Node(5)),
            (Node(5), Node(1)),
        ];
        let g = IndexedDiGraph::from(v);
        let weights = vec![1.0, 3.0, 1.0, 1.0, 2.0, 1.0, 1.0];
        let cycle = g.minimum_mean_weight_cycle(&Node(0), &weights).unwrap();
        assert_eq!(cycle, vec![Edge(3), Edge(4), Edge(5), Edge(6)]);
        assert_eq!(
            g.cycle_as_node_list(&cycle),
            vec![Node(1), Node(3), Node(4), Node(5)]
        );
    }

    /// Example with no loops
    /// 0 ---> 1 -----------+
    /// |                   |
    /// +--> 2 -------> 3 <-+
    ///      |
    ///      +--> 4
    /// There should be no cycles
    #[test]
    fn min_mean_weight_cycle_03() {
        let v = vec![
            (Node(0), Node(1)),
            (Node(0), Node(2)),
            (Node(1), Node(3)),
            (Node(2), Node(3)),
            (Node(2), Node(4)),
        ];
        let g = IndexedDiGraph::from(v);
        let weights = vec![1.0, 1.0, 2.0, 1.0, 2.0];
        let cycle = g.minimum_mean_weight_cycle(&Node(0), &weights);
        assert!(cycle.is_none());
    }

    /// example graph with parallel edges
    ///     +----->
    /// +-> 0 ----> 1 -----> 3
    /// |           |        |
    /// +----- 2 <--+        |
    ///        ^             |
    ///        +-------------+
    #[test]
    fn min_mean_weight_cycle_04_parallel_edges() {
        let v = vec![
            (Node(0), Node(1)),
            (Node(0), Node(1)),
            (Node(1), Node(2)),
            (Node(2), Node(0)),
            (Node(1), Node(3)),
            (Node(3), Node(2)),
        ];
        let g = IndexedDiGraph::from(v);

        let weights = vec![3.0, 1.0, 1.0, 1.0, 1.0, 4.0];
        let cycle = g.minimum_mean_weight_cycle(&Node(0), &weights).unwrap();
        assert_eq!(cycle, vec![Edge(2), Edge(3), Edge(1)]);
        assert_eq!(
            g.cycle_as_node_list(&cycle),
            vec![Node(1), Node(2), Node(0)]
        );
    }

    /// +---------+
    /// |         V
    /// 0 <------ 1
    /// ^         |
    /// +--- 2 <--+
    #[test]
    fn min_mean_weight_cycle_05_small_loop() {
        let v = vec![
            (Node(0), Node(1)),
            (Node(1), Node(0)),
            (Node(1), Node(2)),
            (Node(2), Node(0)),
        ];
        let g = IndexedDiGraph::from(v);

        let weights = vec![1.0, 1.0, 2.0, 3.0];
        let cycle = g.minimum_mean_weight_cycle(&Node(0), &weights).unwrap();
        assert_eq!(cycle, vec![Edge(1), Edge(0)]);
        assert_eq!(g.cycle_as_node_list(&cycle), vec![Node(1), Node(0)]);

        let weights = vec![1.0, 1000.0, 2.0, 3.0];
        let cycle = g.minimum_mean_weight_cycle(&Node(0), &weights).unwrap();
        assert_eq!(cycle, vec![Edge(0), Edge(2), Edge(3)]);
        assert_eq!(
            g.cycle_as_node_list(&cycle),
            vec![Node(0), Node(1), Node(2)]
        );
    }

    /// graph with self loop
    #[test]
    fn min_mean_weight_cycle_06_self_loop() {
        let v = vec![
            (Node(0), Node(1)),
            (Node(1), Node(2)),
            (Node(2), Node(0)),
            (Node(1), Node(1)),
        ];
        let g = IndexedDiGraph::from(v);

        // weight A
        let weights = vec![1.0, 1.0, 1.0, 0.1];
        let cycle = g.minimum_mean_weight_cycle(&Node(0), &weights).unwrap();
        assert_eq!(cycle, vec![Edge(3)]);
        assert_eq!(g.cycle_as_node_list(&cycle), vec![Node(1)]);

        // weight B
        let weights = vec![-1.0, -1.0, -1.0, 0.1];
        let cycle = g.minimum_mean_weight_cycle(&Node(0), &weights).unwrap();
        assert_eq!(cycle, vec![Edge(0), Edge(1), Edge(2)]);
        assert_eq!(
            g.cycle_as_node_list(&cycle),
            vec![Node(0), Node(1), Node(2)]
        );
    }

    #[test]
    fn layout_output() {
        let v = vec![
            (Node(0), Node(1)),
            (Node(0), Node(2)),
            (Node(1), Node(3)),
            (Node(2), Node(3)),
            (Node(2), Node(4)),
        ];
        let g = IndexedDiGraph::from(v);
        let positions = g.layout_by_force_atlas2();
        println!("{:?}", positions)
    }
}
