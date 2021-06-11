use crate::prob::Prob;

struct PHMMParams {
    p_mismatch: Prob,
    p_gap_open: Prob,
    p_gap_ext: Prob,
    p_mm: Prob,
    p_im: Prob,
    n_max_gaps: u32,
}

struct PHMMLayer {
    FM: Vec<Prob>,
    FI: Vec<Prob>,
    FD: Vec<Prob>,
}

#[derive(Debug, PartialEq, PartialOrd, Eq, Hash)]
pub struct Node(usize);

pub trait PHMM {
    // graph structures
    fn nodes(&self) -> Vec<Node>;
    fn childs(&self, v: &Node) -> Vec<Node>;
    fn parents(&self, v: &Node) -> Vec<Node>;
    fn is_adjacent(&self, v: &Node, w: &Node) -> bool;
    fn copy_num(&self, v: &Node) -> u32;
    /*
     */
    /*
    // probs
    fn init_prob(&self, v: &Node) -> Prob;
    fn emission_prob(&self, v: &Node) -> Prob;
    fn trans_prob(&self, v: &Node, w: &Node) -> Prob;

    // prob calculation
    fn forward(&self, emissions: &[u8]) -> Vec<PHMMLayer> {
        let mut layers = Vec::new();
        let layer = self.init();
        layers.push(layer);
        for (i, &emission) in emissions.iter().enumerate() {
            let layer = self.step(&layers[i], emission);
            layers.push(layer);
        }
    }
    fn init() -> PHMMLayer {}
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
        if 0 <= v.0 && v.0 < self.bases.len() {
            1
        } else {
            0
        }
    }
}

struct DbgPHMM {}

// iterator implementation test
struct Counter {
    count: u32,
}
impl Counter {
    fn new() -> Counter {
        Counter { count: 0 }
    }
}
impl Iterator for Counter {
    type Item = u32;
    fn next(&mut self) -> Option<Self::Item> {
        self.count += 1;
        if self.count < 6 {
            Some(self.count)
        } else {
            None
        }
    }
}

fn test2(vec: &[u8]) -> Vec<u8> {
    vec.iter().map(|x| x + 1).collect()
}

pub fn test() {
    println!("hmm test");
    let mut counter = Counter::new();
    for x in counter {
        println!("{:?}", x);
    }

    let vec: Vec<u8> = vec![5, 8, 9];
    let ids: Vec<usize> = (0..vec.len()).collect();
    println!("{:?}, {:?}", vec, ids);

    let model = LinearPHMM::from(b"ATCGATTCGA");
    for node in model.nodes().iter() {
        println!(
            "{:?} {:?} {:?} {:?}",
            node,
            model.childs(node),
            model.parents(node),
            model.copy_num(node),
        );
    }
    println!("{}", model.is_adjacent(&Node(5), &Node(6)));
    println!("{}", model.is_adjacent(&Node(5), &Node(8)));
}
