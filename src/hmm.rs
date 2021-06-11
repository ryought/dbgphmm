use crate::prob::Prob;

#[derive(Debug)]
pub struct PHMMParams {
    p_mismatch: Prob,
    p_gap_open: Prob,
    p_gap_ext: Prob,
    p_MM: Prob,
    p_IM: Prob,
    p_DM: Prob,
    p_MI: Prob,
    p_II: Prob,
    p_DI: Prob,
    p_MD: Prob,
    p_ID: Prob,
    p_DD: Prob,
    n_max_gaps: u32,
}

impl PHMMParams {
    pub fn new(p_mismatch: Prob, p_gap_open: Prob, p_gap_ext: Prob, n_max_gaps: u32) -> PHMMParams {
        PHMMParams {
            p_mismatch,
            p_gap_open,
            p_gap_ext,
            p_DD: p_gap_ext,
            p_II: p_gap_ext,
            p_MI: p_gap_open,
            p_MD: p_gap_open,
            p_ID: p_gap_open,
            p_DI: p_gap_open,
            // p_MM: Prob::from_prob(1.0) - p_gap_open - p_gap_open,
            // p_DM: Prob::from_prob(1.0) - p_gap_open - p_gap_ext,
            p_MM: Prob::from_prob(1.0),
            p_DM: Prob::from_prob(1.0),
            p_IM: Prob::from_prob(1.0),
            n_max_gaps,
        }
    }
}

#[derive(Debug)]
pub struct PHMMLayer {
    FM: Vec<Prob>,
    FI: Vec<Prob>,
    FD: Vec<Prob>,
}

#[derive(Debug, PartialEq, PartialOrd, Eq, Hash)]
pub struct Node(usize);

pub trait PHMM {
    // graph structures
    fn nodes(&self) -> Vec<Node>;
    fn n_nodes(&self) -> usize {
        self.nodes().len()
    }
    fn childs(&self, v: &Node) -> Vec<Node>;
    fn parents(&self, v: &Node) -> Vec<Node>;
    fn is_adjacent(&self, v: &Node, w: &Node) -> bool;
    // hmm related values. each node has (copy_num, emission) attributes
    fn copy_num(&self, v: &Node) -> u32;
    fn total_copy_num(&self) -> u32 {
        self.nodes().iter().map(|v| self.copy_num(v)).sum()
    }
    fn emission(&self, v: &Node) -> u8;
    // TODO this can be computed only from copy_num
    fn trans_prob(&self, v: &Node, w: &Node) -> Prob;
    fn label(&self, v: &Node) -> String {
        let label = String::new();
        label
    }

    fn init_prob(&self, v: &Node) -> Prob {
        let copy_num = self.copy_num(v);
        let total_copy_num = self.total_copy_num();
        Prob::from_prob(f64::from(copy_num) / f64::from(total_copy_num))
    }
    // prob calculation
    fn init(&self, param: &PHMMParams) -> PHMMLayer {
        let mut FM: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
        let mut FI: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
        // (1) FM[i-1], FI[i-1], FD[i-1] -> FM[i], FI[i]
        for v in self.nodes().iter() {
            FM[v.0] = self.init_prob(v);
        }
        // (2) FM[i], FI[i] -> FD[i]
        let mut FDs: Vec<Vec<Prob>> = Vec::new();
        // 0
        let mut FD0: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
        for v in self.nodes().iter() {
            FD0[v.0] = self
                .parents(v)
                .iter()
                .map(|w| self.trans_prob(w, v) * (param.p_MD * FM[w.0] + param.p_ID * FI[w.0]))
                .sum();
        }
        FDs.push(FD0);
        // >0
        for x in 1..param.n_max_gaps {
            let mut FD: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
            let FD_prev = FDs.last().unwrap();
            for v in self.nodes().iter() {
                FD[v.0] = self
                    .parents(v)
                    .iter()
                    .map(|w| self.trans_prob(w, v) * (param.p_DD * FD_prev[w.0]))
                    .sum();
            }
            FDs.push(FD);
        }
        let FD: Vec<Prob> = self
            .nodes()
            .iter()
            .map(|v| FDs.iter().map(|FD| FD[v.0]).sum())
            .collect();
        PHMMLayer { FM, FI, FD }
    }

    /*
    fn forward(&self, emissions: &[u8]) -> Vec<PHMMLayer> {
        let mut layers = Vec::new();
        let layer = self.init();
        layers.push(layer);
        for (i, &emission) in emissions.iter().enumerate() {
            let layer = self.step(&layers[i], emission);
            layers.push(layer);
        }
    }
    fn step(layer_prev: &PHMMLayer, emission: u8) -> PHMMLayer {
        // fill for emission (seq[i])
        for v in self.nodes().iter() {
            // for w in node.childs().iter()
            // or
            for v in self.childs(node).iter() {}
        }
    }
    */
}

pub struct LinearPHMM {
    pub bases: Vec<u8>,
}
impl LinearPHMM {
    pub fn from(seq: &[u8]) -> LinearPHMM {
        LinearPHMM {
            bases: seq.to_vec(),
        }
    }
}
impl PHMM for LinearPHMM {
    fn nodes(&self) -> Vec<Node> {
        (0..self.bases.len()).map(|i| Node(i)).collect()
    }
    fn childs(&self, v: &Node) -> Vec<Node> {
        if v.0 != self.nodes().len() - 1 {
            vec![Node(v.0 + 1)]
        } else {
            vec![]
        }
    }
    fn parents(&self, v: &Node) -> Vec<Node> {
        if v.0 != 0 {
            vec![Node(v.0 - 1)]
        } else {
            vec![]
        }
    }
    fn is_adjacent(&self, v: &Node, w: &Node) -> bool {
        v.0 + 1 == w.0
    }
    fn copy_num(&self, v: &Node) -> u32 {
        if v.0 < self.bases.len() {
            1
        } else {
            0
        }
    }
    fn emission(&self, v: &Node) -> u8 {
        self.bases[v.0]
    }
    fn trans_prob(&self, v: &Node, w: &Node) -> Prob {
        if self.is_adjacent(v, w) {
            Prob::from_prob(1.0)
        } else {
            Prob::from_prob(0.0)
        }
    }
}

struct DbgPHMM {}

fn print(probs: &[Prob]) {
    for (i, p) in probs.iter().enumerate() {
        println!("{}: {}", i, p);
    }
}

pub fn test() {
    println!("hmm test");
    let model = LinearPHMM::from(b"ATCGATTCGA");
    for node in model.nodes().iter() {
        println!(
            "{:?} {:?} {:?} {:?} {:?} {}",
            node,
            model.childs(node),
            model.parents(node),
            model.copy_num(node),
            model.emission(node) as char,
            model.init_prob(node),
        );
    }
    println!(
        "{} {}",
        model.is_adjacent(&Node(5), &Node(6)),
        model.trans_prob(&Node(5), &Node(6))
    );
    println!("{}", model.is_adjacent(&Node(5), &Node(8)));
    let param = PHMMParams::new(
        Prob::from_prob(0.01),
        Prob::from_prob(0.01),
        Prob::from_prob(0.01),
        10,
    );
    println!("param: {:?}", param);
    let l0 = model.init(&param);
    for v in model.nodes().iter() {
        println!(
            "{:?} FM={} FI={} FD={}",
            v, l0.FM[v.0], l0.FI[v.0], l0.FD[v.0]
        );
    }
}
