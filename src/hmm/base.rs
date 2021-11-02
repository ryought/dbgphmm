#![allow(non_snake_case)]
use super::layer::PHMMLayer as PHMMLayerRaw;
use super::params::PHMMParams;
pub use crate::graph::Node;
use crate::prob::Prob;
use crate::veclike::{DenseVec, SparseVec, VecLike};
use arrayvec::ArrayVec;
use itertools::Itertools;
use log::{debug, info, trace, warn};
use std::fmt::Write as FmtWrite;

pub type PHMMLayer = PHMMLayerRaw<DenseVec<Prob>>;

pub fn iter_nodes(n_nodes: usize) -> impl std::iter::Iterator<Item = Node> {
    (0..n_nodes).map(|i| Node(i))
}

/*
// TODO
pub fn iter_active_nodes(
    n_nodes: usize,
    active_nodes: &Option<Vec<Node>>,
) -> Box<dyn std::iter::Iterator<Item = Node>> {
    if let Some(nodes) = active_nodes {
        Box::new(nodes.into_iter())
    } else {
        Box::new(iter_nodes(n_nodes))
    }
}
*/

pub trait PHMM {
    /// return a number of nodes
    fn n_nodes(&self) -> usize;
    fn childs(&self, v: &Node) -> &[Node];
    fn parents(&self, v: &Node) -> &[Node];
    fn is_adjacent(&self, v: &Node, w: &Node) -> bool;
    // hmm related values. each node has (copy_num, emission) attributes
    fn copy_num(&self, v: &Node) -> f64;
    fn total_copy_num(&self) -> f64 {
        iter_nodes(self.n_nodes())
            .filter(|v| self.is_emitable(&v))
            .map(|v| self.copy_num(&v))
            .sum()
    }
    fn emission(&self, v: &Node) -> u8;
    fn is_emitable(&self, v: &Node) -> bool {
        self.emission(v) != b'N'
    }
    fn trans_prob(&self, v: &Node, w: &Node) -> Prob;
    fn label(&self, v: &Node) -> String {
        let label = String::new();
        label
    }
    fn init_prob(&self, v: &Node) -> Prob {
        if self.is_emitable(v) {
            let copy_num = self.copy_num(v);
            let total_copy_num = self.total_copy_num();
            if total_copy_num > 0.0 {
                return Prob::from_prob(copy_num / total_copy_num);
            }
        }
        Prob::from_prob(0.0)
    }
    // forward prob
    fn fmi_init<V: VecLike<Prob>>(&self) -> (V, V) {
        let fM = V::new(self.n_nodes(), Prob::from_prob(0.0));
        let fI = V::new(self.n_nodes(), Prob::from_prob(0.0));
        (fM, fI)
    }
    /// calc (FM[i-1], FI[i-1], FD[i-1]) -> FM[i], FI[i]
    fn fmi_from_fmid<V: VecLike<Prob>>(
        &self,
        param: &PHMMParams,
        layer_t_1: &PHMMLayerRaw<V>,
        layer_t: &mut PHMMLayerRaw<V>,
        emission: u8,
    ) {
        let iterator: Box<dyn std::iter::Iterator<Item = Node>> =
            if let Some(active_nodes) = &layer_t.active_nodes {
                Box::new(active_nodes.iter().map(|&v| v))
            } else {
                Box::new(iter_nodes(self.n_nodes()))
            };
        for v in iterator {
            // M state
            let emission_prob_M: Prob = if self.emission(&v) == emission {
                param.p_match
            } else {
                param.p_mismatch
            };
            let from_normal: Prob = self
                .parents(&v)
                .iter()
                .map(|w| {
                    self.trans_prob(w, &v)
                        * (param.p_MM * layer_t_1.pM.get(w.0)
                            + param.p_IM * layer_t_1.pI.get(w.0)
                            + param.p_DM * layer_t_1.pD.get(w.0))
                })
                .sum();
            let from_begin: Prob =
                self.init_prob(&v) * (param.p_MM * layer_t_1.pMB + param.p_IM * layer_t_1.pIB);
            layer_t
                .pM
                .set(v.0, emission_prob_M * (from_normal + from_begin));

            // I state
            let emission_prob_I: Prob = param.p_random;
            let from_normal: Prob = param.p_MI * layer_t_1.pM.get(v.0)
                + param.p_II * layer_t_1.pI.get(v.0)
                + param.p_DI * layer_t_1.pD.get(v.0);
            layer_t.pI.set(v.0, emission_prob_I * from_normal);
        }
    }
    fn fd_from_fmi<V: VecLike<Prob>>(&self, param: &PHMMParams, layer_t: &mut PHMMLayerRaw<V>) {
        // calc (FM[i], FI[i]) -> FD[i]
        let mut fDs: Vec<V> = Vec::new();
        // 0
        let mut fD0 = V::new(self.n_nodes(), Prob::from_prob(0.0));
        let iterator: Box<dyn std::iter::Iterator<Item = Node>> =
            if let Some(active_nodes) = &layer_t.active_nodes {
                Box::new(active_nodes.iter().map(|&v| v))
            } else {
                Box::new(iter_nodes(self.n_nodes()))
            };
        for v in iterator {
            let from_normal: Prob = self
                .parents(&v)
                .iter()
                .map(|w| {
                    self.trans_prob(w, &v)
                        * (param.p_MD * layer_t.pM.get(w.0) + param.p_ID * layer_t.pI.get(w.0))
                })
                .sum();
            // TODO correct?
            let from_begin: Prob =
                self.init_prob(&v) * (param.p_MD * layer_t.pMB + param.p_ID * layer_t.pIB);
            fD0.set(v.0, from_normal + from_begin);
        }
        fDs.push(fD0);
        // >0
        for x in 1..param.n_max_gaps {
            let mut fD = V::new(self.n_nodes(), Prob::from_prob(0.0));
            let fD_prev = fDs.last().unwrap();
            let iterator: Box<dyn std::iter::Iterator<Item = Node>> =
                if let Some(active_nodes) = &layer_t.active_nodes {
                    Box::new(active_nodes.iter().map(|&v| v))
                } else {
                    Box::new(iter_nodes(self.n_nodes()))
                };
            for v in iterator {
                fD.set(
                    v.0,
                    self.parents(&v)
                        .iter()
                        .map(|w| self.trans_prob(w, &v) * (param.p_DD * fD_prev.get(w.0)))
                        .sum(),
                );
            }
            fDs.push(fD);
        }

        // collect
        let iterator: Box<dyn std::iter::Iterator<Item = Node>> =
            if let Some(active_nodes) = &layer_t.active_nodes {
                Box::new(active_nodes.iter().map(|&v| v))
            } else {
                Box::new(iter_nodes(self.n_nodes()))
            };
        for v in iterator {
            layer_t.pD.set(v.0, fDs.iter().map(|fD| fD.get(v.0)).sum());
        }
    }
    fn fb_init(&self) -> (Prob, Prob) {
        let fMB = Prob::from_prob(1.0);
        let fIB = Prob::from_prob(0.0);
        (fMB, fIB)
    }
    fn fb_from_fb<V: VecLike<Prob>>(
        &self,
        param: &PHMMParams,
        layer_t_1: &PHMMLayerRaw<V>,
        layer_t: &mut PHMMLayerRaw<V>,
    ) {
        layer_t.pMB = Prob::from_prob(0.0);
        let emission_prob: Prob = Prob::from_prob(0.25);
        layer_t.pIB = emission_prob * (param.p_MI * layer_t_1.pMB + param.p_II * layer_t_1.pIB);
    }
    fn fe_init(&self) -> Prob {
        let fE = Prob::from_prob(0.0);
        fE
    }
    fn fe_from_fmid<V: VecLike<Prob>>(&self, param: &PHMMParams, layer_t: &mut PHMMLayerRaw<V>) {
        let iterator: Box<dyn std::iter::Iterator<Item = Node>> =
            if let Some(active_nodes) = &layer_t.active_nodes {
                Box::new(active_nodes.iter().map(|&v| v))
            } else {
                Box::new(iter_nodes(self.n_nodes()))
            };
        layer_t.pE = iterator
            .map(|v| {
                param.p_end * (layer_t.pM.get(v.0) + layer_t.pI.get(v.0) + layer_t.pD.get(v.0))
            })
            .sum();
    }
    // prob calculation
    fn f_init<V: VecLike<Prob>>(&self, param: &PHMMParams) -> PHMMLayerRaw<V> {
        let mut layer = PHMMLayerRaw::f_init(self.n_nodes());
        self.fd_from_fmi(param, &mut layer);
        layer
    }
    /// compute layer[t].active_nodes from layer[t-1].active_nodes
    /// by taking active node's childs
    fn get_next_active_nodes(&self, active_nodes: &[Node]) -> Vec<Node> {
        active_nodes
            .iter()
            .flat_map(|v| self.childs(v))
            .map(|&w| w)
            .unique()
            .collect()
    }
    /// given a list of previously active nodes,
    /// calculates the all candidates of next active nodes
    /// that is all childrens of prev active nodes
    fn upconvert_active_nodes<V: VecLike<Prob>>(
        &self,
        prev_layer: &PHMMLayerRaw<V>,
    ) -> Option<Vec<Node>> {
        if let Some(nodes) = &prev_layer.active_nodes {
            let next_active_nodes = self.get_next_active_nodes(nodes);
            let labels: Vec<String> = next_active_nodes.iter().map(|v| self.label(v)).collect();
            info!("upconvert={:?}", labels);
            Some(next_active_nodes)
        } else {
            info!("upconvert=None");
            None
        }
    }
    /// given a current PHMMLayer and its active nodes,
    /// compute a list of the prominent active nodes
    fn refine_active_nodes<V: VecLike<Prob>>(
        &self,
        layer: &PHMMLayerRaw<V>,
        max_active_nodes: usize,
    ) -> Option<Vec<Node>> {
        // calculate layer kmer probabilities
        let mut scores: Vec<(usize, f64)> = layer
            .to_kmer_prob()
            .into_iter()
            .map(|p| p.to_log_value())
            .enumerate()
            .collect();
        // sort descending(from largest to smallest) order
        scores.sort_by(|(_, a), (_, b)| b.partial_cmp(a).unwrap());
        for (i, p) in scores.iter().take(max_active_nodes * 2) {
            info!("refine i={} label={} p={}", i, self.label(&Node(*i)), p);
        }
        Some(
            scores
                .iter()
                .take(max_active_nodes)
                .map(|(i, _)| Node(*i))
                .collect(),
        )
    }
    fn f_step<V: VecLike<Prob>>(
        &self,
        param: &PHMMParams,
        prev_layer: &PHMMLayerRaw<V>,
        emission: u8,
        index: usize,
    ) -> PHMMLayerRaw<V> {
        let mut next_layer = PHMMLayerRaw::new(self.n_nodes());
        // compute candidate active_nodes from prev_layer and store it in next_layer
        if param.only_active_nodes && index > param.n_ignore_active_nodes_first {
            next_layer.active_nodes = self.upconvert_active_nodes(&prev_layer);
            info!(
                "f_step({}) upconverted {:?}",
                index, next_layer.active_nodes
            );
        }
        self.fmi_from_fmid(param, &prev_layer, &mut next_layer, emission);
        self.fb_from_fb(param, &prev_layer, &mut next_layer);
        self.fd_from_fmi(param, &mut next_layer);
        self.fe_from_fmid(param, &mut next_layer);
        // reduce active nodes
        if param.only_active_nodes && index > param.n_ignore_active_nodes_first {
            next_layer.active_nodes =
                self.refine_active_nodes(&next_layer, param.n_max_active_nodes);
            info!("f_step({}) refined {:?}", index, next_layer.active_nodes);
        }
        next_layer
    }
    ///
    /// Forward probability
    /// layers[i] = P(x[0..i], end with a state)
    ///
    fn forward(&self, param: &PHMMParams, emissions: &[u8]) -> Vec<PHMMLayer> {
        trace!("start forward!");
        let mut layers = Vec::new();
        let layer = self.f_init::<DenseVec<Prob>>(&param);
        info!("l0:{}", layer.pM.len());
        trace!("l0:\n{}", layer);
        layers.push(layer);
        for (i, &emission) in emissions.iter().enumerate() {
            trace!("l{} emission={}", i, emission as char);
            let layer = self.f_step(&param, &layers[i], emission, i);
            trace!("l{}:\n{}", i, layer);
            info!("pE({})={}", i, layer.pE);
            layers.push(layer);
        }
        layers
    }
    fn forward_sparse(
        &self,
        param: &PHMMParams,
        emissions: &[u8],
    ) -> Vec<PHMMLayerRaw<SparseVec<Prob>>> {
        let mut layers = Vec::new();
        let layer = self.f_init::<SparseVec<Prob>>(&param);
        info!("l0:{}", layer.pM.len());
        layers.push(layer);
        for (i, &emission) in emissions.iter().enumerate() {
            let layer = self.f_step(&param, &layers[i], emission, i);
            info!("pE({})={}", i, layer.pE);
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
    fn bd_from_bmi<V: VecLike<Prob>>(
        &self,
        param: &PHMMParams,
        bm1: &V,
        bi1: &V,
        emission: u8,
    ) -> DenseVec<Prob> {
        let mut bDs: Vec<V> = Vec::new();
        // 0
        let mut bD0 = V::new(self.n_nodes(), Prob::from_prob(0.0));
        for v in iter_nodes(self.n_nodes()) {
            let to_m: Prob = self
                .childs(&v)
                .iter()
                .map(|w| {
                    let emission_prob_M: Prob = if self.emission(w) == emission {
                        param.p_match
                    } else {
                        param.p_mismatch
                    };
                    self.trans_prob(&v, w) * param.p_DM * bm1.get(w.0) * emission_prob_M
                })
                .sum();
            let emission_prob_I: Prob = param.p_random;
            let to_i = param.p_DI * emission_prob_I * bi1.get(v.0);
            bD0.set(v.0, to_m + to_i);
        }
        bDs.push(bD0);
        // >0
        for x in 1..param.n_max_gaps {
            let mut bD = V::new(self.n_nodes(), Prob::from_prob(0.0));
            let bD_prev = bDs.last().unwrap();
            for v in iter_nodes(self.n_nodes()) {
                bD.set(
                    v.0,
                    self.childs(&v)
                        .iter()
                        .map(|w| self.trans_prob(&v, w) * param.p_DD * bD_prev.get(w.0))
                        .sum(),
                );
            }
            bDs.push(bD);
        }
        let bD = iter_nodes(self.n_nodes())
            .map(|v| bDs.iter().map(|bD| bD.get(v.0)).sum())
            .collect();
        bD
    }
    fn bmi_from_bmid<V: VecLike<Prob>>(
        &self,
        param: &PHMMParams,
        bm1: &V,
        bi1: &V,
        bd0: &V,
        emission: u8,
    ) -> (V, V) {
        let mut bm0 = V::new(self.n_nodes(), Prob::from_prob(0.0));
        let mut bi0 = V::new(self.n_nodes(), Prob::from_prob(0.0));
        for v in iter_nodes(self.n_nodes()) {
            // fill bm0
            let to_md: Prob = self
                .childs(&v)
                .iter()
                .map(|w| {
                    let emission_prob_M: Prob = if self.emission(w) == emission {
                        param.p_match
                    } else {
                        param.p_mismatch
                    };
                    self.trans_prob(&v, w)
                        * (param.p_MM * emission_prob_M * bm1.get(w.0) + param.p_MD * bd0.get(w.0))
                })
                .sum();
            let emission_prob_I: Prob = param.p_random;
            let to_i = param.p_MI * emission_prob_I * bi1.get(v.0);
            bm0.set(v.0, to_md + to_i);

            // fill bi0
            let to_md: Prob = self
                .childs(&v)
                .iter()
                .map(|w| {
                    let emission_prob_M: Prob = if self.emission(w) == emission {
                        param.p_match
                    } else {
                        param.p_mismatch
                    };
                    self.trans_prob(&v, w)
                        * (param.p_IM * emission_prob_M * bm1.get(w.0) + param.p_ID * bd0.get(w.0))
                })
                .sum();
            let emission_prob_I: Prob = param.p_random;
            let to_i = param.p_II * emission_prob_I * bi1.get(v.0);
            bi0.set(v.0, to_md + to_i);
        }
        (bm0, bi0)
    }
    fn bb_from_bmd<V: VecLike<Prob>>(
        &self,
        param: &PHMMParams,
        bm1: &V,
        bd0: &V,
        bib1: Prob,
        emission: u8,
    ) -> (Prob, Prob) {
        // bmb0
        let to_md: Prob = iter_nodes(self.n_nodes())
            .map(|v| {
                let emission_prob_M: Prob = if self.emission(&v) == emission {
                    param.p_match
                } else {
                    param.p_mismatch
                };
                self.init_prob(&v)
                    * (param.p_MM * emission_prob_M * bm1.get(v.0) + param.p_MD * bd0.get(v.0))
            })
            .sum();
        let emission_prob_I: Prob = param.p_random;
        let to_i = param.p_MI * emission_prob_I * bib1;
        let bmb0 = to_md + to_i;

        // bib0
        let to_md: Prob = iter_nodes(self.n_nodes())
            .map(|v| {
                let emission_prob_M: Prob = if self.emission(&v) == emission {
                    param.p_match
                } else {
                    param.p_mismatch
                };
                self.init_prob(&v)
                    * (param.p_IM * emission_prob_M * bm1.get(v.0) + param.p_ID * bd0.get(v.0))
            })
            .sum();
        let emission_prob_I: Prob = param.p_random;
        let to_i = param.p_II * emission_prob_I * bib1;
        let bib0 = to_md + to_i;

        (bmb0, bib0)
    }
    fn b_init(&self, param: &PHMMParams) -> PHMMLayer {
        let mut bM = DenseVec::new(self.n_nodes(), param.p_end);
        let mut bI = DenseVec::new(self.n_nodes(), param.p_end);
        let mut bD = DenseVec::new(self.n_nodes(), param.p_end);
        PHMMLayer {
            pM: bM,
            pI: bI,
            pD: bD,
            pMB: Prob::from_prob(0.0),
            pIB: Prob::from_prob(0.0),
            pE: Prob::from_prob(0.0),
            active_nodes: None,
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
            active_nodes: None,
        }
    }
    /// Backward probability
    /// layers[i] = P(x[i..] | start with a state)
    fn backward(&self, param: &PHMMParams, emissions: &[u8]) -> Vec<PHMMLayer> {
        debug!("start backward!");
        let mut layers = Vec::new();
        let layer = self.b_init(&param);
        trace!("l0:{}", layer.pM.len());
        trace!("l0:\n{}", layer);
        layers.push(layer);
        for (i, &emission) in emissions.iter().rev().enumerate() {
            trace!("l{} emission={}", i, emission as char);
            let layer = self.b_step(&param, &layers[i], emission);
            trace!("l{}:\n{}", i, layer);
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
    /// p[i][j] = (e[i] is emitted from j-th state)
    fn state_prob(
        &self,
        forward_layers: &[PHMMLayer],
        backward_layers: &[PHMMLayer],
    ) -> Vec<PHMMLayer> {
        let mut prob_layers = Vec::new();
        let p = backward_layers[0].pMB;
        for (fl, bl) in forward_layers.iter().zip(backward_layers.iter()) {
            let pl = PHMMLayer {
                pM: fl
                    .pM
                    .iter()
                    .zip(bl.pM.iter())
                    .map(|(f, b)| f * b / p)
                    .collect(),
                pI: fl
                    .pI
                    .iter()
                    .zip(bl.pI.iter())
                    .map(|(f, b)| f * b / p)
                    .collect(),
                pD: fl
                    .pD
                    .iter()
                    .zip(bl.pD.iter())
                    .map(|(f, b)| f * b / p)
                    .collect(),
                pMB: fl.pMB * bl.pMB / p,
                pIB: fl.pIB * bl.pIB / p,
                pE: fl.pE * bl.pE / p,
                active_nodes: None,
            };
            prob_layers.push(pl);
        }
        prob_layers
    }
    // output
    fn as_dot(&self) -> String {
        let mut s = String::new();
        writeln!(&mut s, "digraph dbgphmm {{");
        for v in iter_nodes(self.n_nodes()) {
            // for node
            let emission = self.emission(&v);
            let copy_num = self.copy_num(&v);
            let init_prob = self.init_prob(&v);
            writeln!(
                &mut s,
                "\t{} [label=\"{} x{}\"];",
                v.0,
                emission as char,
                copy_num // "\t{} [label=\"{} x{} {}\"];",
                         // v.0, emission as char, copy_num, init_prob
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
        for v in iter_nodes(self.n_nodes()) {
            writeln!(
                &mut s,
                "{:?} {:?} {:?} {} {:?}",
                v,
                self.childs(&v),
                self.parents(&v),
                self.copy_num(&v),
                self.emission(&v) as char,
            );
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;
}
