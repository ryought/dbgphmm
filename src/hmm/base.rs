use super::params::PHMMParams;
use crate::prob::Prob;
use log::info;
use std::fmt::Write as FmtWrite;

#[derive(Debug)]
pub struct PHMMLayer {
    FM: Vec<Prob>,
    FI: Vec<Prob>,
    FD: Vec<Prob>,
    FMB: Prob,
    FIB: Prob,
}

impl std::fmt::Display for PHMMLayer {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "Node:Begin FM={} FI={}", self.FMB, self.FIB);
        for i in 0..self.FD.len() {
            writeln!(
                f,
                "Node:{}\tFM={}\tFI={}\tFD={}",
                i, self.FM[i], self.FI[i], self.FD[i]
            );
        }
        Ok(())
    }
}

#[derive(Debug, PartialEq, PartialOrd, Eq, Hash, Copy, Clone)]
pub struct Node(pub usize);

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
    // forward prob
    fn fmi_init(&self) -> (Vec<Prob>, Vec<Prob>) {
        let mut FM: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
        let mut FI: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
        (FM, FI)
    }
    fn fmi_from_fmid(
        &self,
        param: &PHMMParams,
        fm: &[Prob],
        fi: &[Prob],
        fd: &[Prob],
        fmb: Prob,
        fib: Prob,
        emission: u8,
    ) -> (Vec<Prob>, Vec<Prob>) {
        ///
        /// calc (FM[i-1], FI[i-1], FD[i-1]) -> FM[i], FI[i]
        ///
        let mut FM: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
        let mut FI: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
        for v in self.nodes().iter() {
            // M state
            let emission_prob_M: Prob = if self.emission(v) == emission {
                param.p_match
            } else {
                param.p_mismatch
            };
            let from_normal: Prob = self
                .parents(v)
                .iter()
                .map(|w| {
                    self.trans_prob(w, v)
                        * (param.p_MM * fm[w.0] + param.p_IM * fi[w.0] + param.p_DM * fd[w.0])
                })
                .sum();
            let from_begin: Prob = self.init_prob(v) * (param.p_MM * fmb + param.p_IM * fib);
            FM[v.0] = emission_prob_M * (from_normal + from_begin);

            // I state
            let emission_prob_I: Prob = param.p_random;
            let from_normal: Prob =
                param.p_MI * fm[v.0] + param.p_II * fi[v.0] + param.p_DI * fd[v.0];
            FI[v.0] = emission_prob_I * from_normal;
        }
        (FM, FI)
    }
    fn fd_from_fmi(
        &self,
        param: &PHMMParams,
        fm: &[Prob],
        fi: &[Prob],
        fmb: Prob,
        fib: Prob,
    ) -> Vec<Prob> {
        ///
        /// calc (FM[i], FI[i]) -> FD[i]
        ///
        let mut FDs: Vec<Vec<Prob>> = Vec::new();
        // 0
        let mut FD0: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
        for v in self.nodes().iter() {
            let from_normal: Prob = self
                .parents(v)
                .iter()
                .map(|w| self.trans_prob(w, v) * (param.p_MD * fm[w.0] + param.p_ID * fi[w.0]))
                .sum();
            let from_begin: Prob = self.init_prob(v) * (param.p_MD * fmb + param.p_ID * fib);
            FD0[v.0] = from_normal + from_begin;
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
        FD
    }
    fn fb_init(&self) -> (Prob, Prob) {
        let FMB = Prob::from_prob(1.0);
        let FIB = Prob::from_prob(0.0);
        (FMB, FIB)
    }
    fn fb_from_fb(&self, param: &PHMMParams, fmb: Prob, fib: Prob) -> (Prob, Prob) {
        let FMB = Prob::from_prob(0.0);
        let emission_prob_I: Prob = Prob::from_prob(0.25);
        let FIB = emission_prob_I * (param.p_MI * fmb + param.p_II * fib);
        (FMB, FIB)
    }
    // prob calculation
    fn init(&self, param: &PHMMParams) -> PHMMLayer {
        let (FM, FI) = self.fmi_init();
        let (FMB, FIB) = self.fb_init();
        let FD = self.fd_from_fmi(param, &FM, &FI, FMB, FIB);
        PHMMLayer {
            FM,
            FI,
            FD,
            FMB,
            FIB,
        }
    }
    fn step(&self, param: &PHMMParams, prev_layer: &PHMMLayer, emission: u8) -> PHMMLayer {
        let (FM, FI) = self.fmi_from_fmid(
            param,
            &prev_layer.FM,
            &prev_layer.FI,
            &prev_layer.FD,
            prev_layer.FMB,
            prev_layer.FIB,
            emission,
        );
        let FD = self.fd_from_fmi(param, &FM, &FI, prev_layer.FMB, prev_layer.FIB);
        let (FMB, FIB) = self.fb_from_fb(param, prev_layer.FMB, prev_layer.FIB);
        PHMMLayer {
            FM,
            FI,
            FD,
            FMB,
            FIB,
        }
    }
    fn forward(&self, param: &PHMMParams, emissions: &[u8]) -> Vec<PHMMLayer> {
        let mut layers = Vec::new();
        let layer = self.init(&param);
        info!("l0:\n{}", layer);
        layers.push(layer);
        for (i, &emission) in emissions.iter().enumerate() {
            let layer = self.step(&param, &layers[i], emission);
            info!("l{}:\n{}", i, layer);
            layers.push(layer);
        }
        layers
    }
    fn forward_prob(&self, param: &PHMMParams, emissions: &[u8]) -> Prob {
        let layers = self.forward(param, emissions);
        let last_layer = layers.last().unwrap();
        let SFM: Prob = last_layer.FM.iter().sum();
        let SFI: Prob = last_layer.FI.iter().sum();
        let SFD: Prob = last_layer.FD.iter().sum();
        SFM + SFI + SFD + last_layer.FMB + last_layer.FIB
    }

    // output
    fn as_dot(&self) -> String {
        let mut s = String::new();
        writeln!(&mut s, "digraph dbgphmm {{");
        for v in self.nodes().iter() {
            // for node
            let emission = self.emission(&v);
            let copy_num = self.copy_num(&v);
            writeln!(
                &mut s,
                "\t{} [label=\"{} x{}\"];",
                v.0, emission as char, copy_num
            );
            // for edges
            for w in self.childs(&v).iter() {
                let p = self.trans_prob(&v, &w);
                writeln!(&mut s, "\t{} -> {} [label=\"{}\"];", v.0, w.0, p);
            }
        }
        writeln!(&mut s, "}}");
        s
    }
}
