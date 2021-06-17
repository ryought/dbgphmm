use super::params::PHMMParams;
pub use crate::graph::Node;
use crate::prob::Prob;
use arrayvec::ArrayVec;
use log::{debug, info, warn};
use std::fmt::Write as FmtWrite;

#[derive(Debug)]
pub struct PHMMLayer {
    pM: Vec<Prob>,
    pI: Vec<Prob>,
    pD: Vec<Prob>,
    // Begin
    pMB: Prob,
    pIB: Prob,
    // End
    pE: Prob,
}

impl std::fmt::Display for PHMMLayer {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "Node:Begin\tpM={} pI={}", self.pMB, self.pIB);
        for i in 0..self.pD.len() {
            writeln!(
                f,
                "Node:{}\tpM={}\tpI={}\tpD={}",
                i, self.pM[i], self.pI[i], self.pD[i]
            );
        }
        writeln!(f, "Node:End\tpE={}", self.pE);
        Ok(())
    }
}

pub trait PHMM {
    // graph structures
    fn nodes(&self) -> &[Node];
    fn n_nodes(&self) -> usize {
        self.nodes().len()
    }
    fn childs(&self, v: &Node) -> &[Node];
    fn parents(&self, v: &Node) -> &[Node];
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
        let mut fM: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
        let mut fI: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
        (fM, fI)
    }
    /// calc (FM[i-1], FI[i-1], FD[i-1]) -> FM[i], FI[i]
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
        let mut fM: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
        let mut fI: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
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
            fM[v.0] = emission_prob_M * (from_normal + from_begin);

            // I state
            let emission_prob_I: Prob = param.p_random;
            let from_normal: Prob =
                param.p_MI * fm[v.0] + param.p_II * fi[v.0] + param.p_DI * fd[v.0];
            fI[v.0] = emission_prob_I * from_normal;
        }
        (fM, fI)
    }
    fn fd_from_fmi(
        &self,
        param: &PHMMParams,
        fm: &[Prob],
        fi: &[Prob],
        fmb: Prob,
        fib: Prob,
    ) -> Vec<Prob> {
        // calc (FM[i], FI[i]) -> FD[i]
        let mut fDs: Vec<Vec<Prob>> = Vec::new();
        // 0
        let mut fD0: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
        warn!("fd0 init {}", self.nodes().len());
        for v in self.nodes().iter() {
            let from_normal: Prob = self
                .parents(v)
                .iter()
                .map(|w| self.trans_prob(w, v) * (param.p_MD * fm[w.0] + param.p_ID * fi[w.0]))
                .sum();
            let from_begin: Prob = self.init_prob(v) * (param.p_MD * fmb + param.p_ID * fib);
            fD0[v.0] = from_normal + from_begin;
        }
        fDs.push(fD0);
        // >0
        for x in 1..param.n_max_gaps {
            warn!("fd{} init", x);
            let mut fD: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
            let fD_prev = fDs.last().unwrap();
            for v in self.nodes().iter() {
                fD[v.0] = self
                    .parents(v)
                    .iter()
                    .map(|w| self.trans_prob(w, v) * (param.p_DD * fD_prev[w.0]))
                    .sum();
            }
            fDs.push(fD);
        }
        let fD: Vec<Prob> = self
            .nodes()
            .iter()
            .map(|v| fDs.iter().map(|fD| fD[v.0]).sum())
            .collect();
        fD
    }
    fn fb_init(&self) -> (Prob, Prob) {
        let fMB = Prob::from_prob(1.0);
        let fIB = Prob::from_prob(0.0);
        (fMB, fIB)
    }
    fn fb_from_fb(&self, param: &PHMMParams, fmb: Prob, fib: Prob) -> (Prob, Prob) {
        let fMB = Prob::from_prob(0.0);
        let emission_prob_I: Prob = Prob::from_prob(0.25);
        let fIB = emission_prob_I * (param.p_MI * fmb + param.p_II * fib);
        (fMB, fIB)
    }
    fn fe_init(&self) -> Prob {
        let fE = Prob::from_prob(0.0);
        fE
    }
    fn fe_from_fmid(&self, param: &PHMMParams, fm: &[Prob], fi: &[Prob], fd: &[Prob]) -> Prob {
        self.nodes()
            .iter()
            .map(|v| param.p_end * (fm[v.0] + fi[v.0] + fd[v.0]))
            .sum()
    }
    // prob calculation
    fn f_init(&self, param: &PHMMParams) -> PHMMLayer {
        let (fM, fI) = self.fmi_init();
        let (fMB, fIB) = self.fb_init();
        let fE = self.fe_init();
        let fD = self.fd_from_fmi(param, &fM, &fI, fMB, fIB);
        PHMMLayer {
            pM: fM,
            pI: fI,
            pD: fD,
            pMB: fMB,
            pIB: fIB,
            pE: fE,
        }
    }
    fn f_step(&self, param: &PHMMParams, prev_layer: &PHMMLayer, emission: u8) -> PHMMLayer {
        let (fM, fI) = self.fmi_from_fmid(
            param,
            &prev_layer.pM,
            &prev_layer.pI,
            &prev_layer.pD,
            prev_layer.pMB,
            prev_layer.pIB,
            emission,
        );
        let fD = self.fd_from_fmi(param, &fM, &fI, prev_layer.pMB, prev_layer.pIB);
        let (fMB, fIB) = self.fb_from_fb(param, prev_layer.pMB, prev_layer.pIB);
        let fE = self.fe_from_fmid(param, &fM, &fI, &fD);
        PHMMLayer {
            pM: fM,
            pI: fI,
            pD: fD,
            pMB: fMB,
            pIB: fIB,
            pE: fE,
        }
    }
    fn forward(&self, param: &PHMMParams, emissions: &[u8]) -> Vec<PHMMLayer> {
        debug!("start forward!");
        let mut layers = Vec::new();
        let layer = self.f_init(&param);
        info!("l0:{}", layer.pM.len());
        debug!("l0:\n{}", layer);
        layers.push(layer);
        for (i, &emission) in emissions.iter().enumerate() {
            info!("l{} emission={}", i, emission as char);
            let layer = self.f_step(&param, &layers[i], emission);
            debug!("l{}:\n{}", i, layer);
            layers.push(layer);
        }
        layers
    }
    fn forward_prob(&self, param: &PHMMParams, emissions: &[u8]) -> Prob {
        let layers = self.forward(param, emissions);
        let last_layer = layers.last().unwrap();
        last_layer.pE
    }

    // backward
    fn bd_from_bmi(
        &self,
        param: &PHMMParams,
        bm1: &[Prob],
        bi1: &[Prob],
        emission: u8,
    ) -> Vec<Prob> {
        let mut bDs: Vec<Vec<Prob>> = Vec::new();
        // 0
        let mut bD0: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
        for v in self.nodes().iter() {
            let to_m: Prob = self
                .childs(v)
                .iter()
                .map(|w| {
                    let emission_prob_M: Prob = if self.emission(w) == emission {
                        param.p_match
                    } else {
                        param.p_mismatch
                    };
                    self.trans_prob(v, w) * param.p_DM * bm1[w.0] * emission_prob_M
                })
                .sum();
            let emission_prob_I: Prob = param.p_random;
            let to_i = param.p_DI * emission_prob_I * bi1[v.0];
            bD0[v.0] = to_m + to_i;
        }
        bDs.push(bD0);
        // >0
        for x in 1..param.n_max_gaps {
            let mut bD: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
            let bD_prev = bDs.last().unwrap();
            for v in self.nodes().iter() {
                bD[v.0] = self
                    .childs(v)
                    .iter()
                    .map(|w| self.trans_prob(v, w) * param.p_DD * bD_prev[w.0])
                    .sum();
            }
            bDs.push(bD);
        }
        let bD: Vec<Prob> = self
            .nodes()
            .iter()
            .map(|v| bDs.iter().map(|bD| bD[v.0]).sum())
            .collect();
        bD
    }
    fn bmi_from_bmid(
        &self,
        param: &PHMMParams,
        bm1: &[Prob],
        bi1: &[Prob],
        bd0: &[Prob],
        emission: u8,
    ) -> (Vec<Prob>, Vec<Prob>) {
        let mut bm0: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
        let mut bi0: Vec<Prob> = vec![Prob::from_prob(0.0); self.n_nodes()];
        for v in self.nodes().iter() {
            // fill bm0
            let to_md: Prob = self
                .childs(v)
                .iter()
                .map(|w| {
                    let emission_prob_M: Prob = if self.emission(w) == emission {
                        param.p_match
                    } else {
                        param.p_mismatch
                    };
                    self.trans_prob(v, w)
                        * (param.p_MM * emission_prob_M * bm1[w.0] + param.p_MD * bd0[w.0])
                })
                .sum();
            let emission_prob_I: Prob = param.p_random;
            let to_i = param.p_MI * emission_prob_I * bi1[v.0];
            bm0[v.0] = to_md + to_i;

            // fill bi0
            let to_md: Prob = self
                .childs(v)
                .iter()
                .map(|w| {
                    let emission_prob_M: Prob = if self.emission(w) == emission {
                        param.p_match
                    } else {
                        param.p_mismatch
                    };
                    self.trans_prob(v, w)
                        * (param.p_IM * emission_prob_M * bm1[w.0] + param.p_ID * bd0[w.0])
                })
                .sum();
            let emission_prob_I: Prob = param.p_random;
            let to_i = param.p_II * emission_prob_I * bi1[v.0];
            bi0[v.0] = to_md + to_i;
        }
        (bm0, bi0)
    }
    fn bb_from_bmd(
        &self,
        param: &PHMMParams,
        bm1: &[Prob],
        bd0: &[Prob],
        bib1: Prob,
        emission: u8,
    ) -> (Prob, Prob) {
        // bmb0
        let to_md: Prob = self
            .nodes()
            .iter()
            .map(|v| {
                let emission_prob_M: Prob = if self.emission(v) == emission {
                    param.p_match
                } else {
                    param.p_mismatch
                };
                self.init_prob(v)
                    * (param.p_MM * emission_prob_M * bm1[v.0] + param.p_MD * bd0[v.0])
            })
            .sum();
        let emission_prob_I: Prob = param.p_random;
        let to_i = param.p_MI * emission_prob_I * bib1;
        let bmb0 = to_md + to_i;

        // bib0
        let to_md: Prob = self
            .nodes()
            .iter()
            .map(|v| {
                let emission_prob_M: Prob = if self.emission(v) == emission {
                    param.p_match
                } else {
                    param.p_mismatch
                };
                self.init_prob(v)
                    * (param.p_IM * emission_prob_M * bm1[v.0] + param.p_ID * bd0[v.0])
            })
            .sum();
        let emission_prob_I: Prob = param.p_random;
        let to_i = param.p_II * emission_prob_I * bib1;
        let bib0 = to_md + to_i;

        (bmb0, bib0)
    }
    fn b_init(&self, param: &PHMMParams) -> PHMMLayer {
        let mut bM: Vec<Prob> = vec![param.p_end; self.n_nodes()];
        let mut bI: Vec<Prob> = vec![param.p_end; self.n_nodes()];
        let mut bD: Vec<Prob> = vec![param.p_end; self.n_nodes()];
        PHMMLayer {
            pM: bM,
            pI: bI,
            pD: bD,
            pMB: Prob::from_prob(0.0),
            pIB: Prob::from_prob(0.0),
            pE: Prob::from_prob(0.0),
        }
    }
    fn b_step(&self, param: &PHMMParams, prev_layer: &PHMMLayer, emission: u8) -> PHMMLayer {
        // emission: x[i+1]
        let bD = self.bd_from_bmi(param, &prev_layer.pM, &prev_layer.pI, emission);
        let (bM, bI) = self.bmi_from_bmid(param, &prev_layer.pM, &prev_layer.pI, &bD, emission);
        let (bMB, bIB) = self.bb_from_bmd(param, &prev_layer.pM, &bD, prev_layer.pIB, emission);
        PHMMLayer {
            pM: bM,
            pI: bI,
            pD: bD,
            pMB: bMB,
            pIB: bIB,
            pE: Prob::from_prob(0.0),
        }
    }
    fn backward(&self, param: &PHMMParams, emissions: &[u8]) -> Vec<PHMMLayer> {
        debug!("start backward!");
        let mut layers = Vec::new();
        let layer = self.b_init(&param);
        info!("l0:{}", layer.pM.len());
        debug!("l0:\n{}", layer);
        layers.push(layer);
        for (i, &emission) in emissions.iter().rev().enumerate() {
            info!("l{} emission={}", i, emission as char);
            let layer = self.b_step(&param, &layers[i], emission);
            debug!("l{}:\n{}", i, layer);
            layers.push(layer);
        }
        layers.reverse();
        layers
    }
    fn backward_prob(&self, param: &PHMMParams, emissions: &[u8]) -> Prob {
        let layers = self.backward(param, emissions);
        let first_layer = layers.first().unwrap();
        first_layer.pMB
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
    fn as_node_list(&self) -> String {
        let mut s = String::new();
        for v in self.nodes().iter() {
            writeln!(
                &mut s,
                "{:?} {:?} {:?} {} {:?}",
                v,
                self.childs(v),
                self.parents(v),
                self.copy_num(v),
                self.emission(v) as char,
            );
        }
        s
    }
}
