#![allow(non_snake_case)]
use super::layer::PHMMLayer as PHMMLayerRaw;
use super::params::PHMMParams;
pub use crate::graph::Node;
use crate::prob::Prob;
use crate::veclike::{DenseVec, SparseVec, VecLike};
use itertools::Itertools;
use log::{debug, info, trace, warn};
use std::fmt::Write as FmtWrite;

pub type PHMMLayer = PHMMLayerRaw<DenseVec<Prob>>;

pub fn iter_nodes(n_nodes: usize) -> impl std::iter::Iterator<Item = Node> {
    (0..n_nodes).map(|i| Node(i))
}

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
    fn label(&self, _: &Node) -> String {
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
    // output
    #[allow(unused_must_use)]
    fn as_dot(&self) -> String {
        let mut s = String::new();
        writeln!(&mut s, "digraph dbgphmm {{");
        for v in iter_nodes(self.n_nodes()) {
            // for node
            let emission = self.emission(&v);
            let copy_num = self.copy_num(&v);
            let _init_prob = self.init_prob(&v);
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
    #[allow(unused_must_use)]
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

pub trait PHMMForward: PHMM {
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
        // 0
        let mut fD0 = V::new(self.n_nodes(), Prob::from_prob(0.0));
        let mut fD1;
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
            let from_begin: Prob =
                self.init_prob(&v) * (param.p_MD * layer_t.pMB + param.p_ID * layer_t.pIB);
            let fD0v = from_normal + from_begin;
            fD0.set(v.0, fD0v);
            layer_t.pD.set(v.0, fD0v);
        }
        // >0
        for _x in 1..param.n_max_gaps {
            // calculate new fD1 from fD0
            fD1 = V::new(self.n_nodes(), Prob::from_prob(0.0));
            let iterator: Box<dyn std::iter::Iterator<Item = Node>> =
                if let Some(active_nodes) = &layer_t.active_nodes {
                    Box::new(active_nodes.iter().map(|&v| v))
                } else {
                    Box::new(iter_nodes(self.n_nodes()))
                };
            for v in iterator {
                let fD1v = self
                    .parents(&v)
                    .iter()
                    .map(|w| self.trans_prob(w, &v) * (param.p_DD * fD0.get(w.0)))
                    .sum();
                fD1.set(v.0, fD1v);
                layer_t.pD.set(v.0, layer_t.pD.get(v.0) + fD1v);
            }
            fD0 = fD1;
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
    fn f_get_next_active_nodes(&self, active_nodes: &[Node]) -> Vec<Node> {
        active_nodes
            .iter()
            .flat_map(|v| self.childs(v))
            .map(|&w| w)
            .unique()
            .collect()
    }
    fn b_get_next_active_nodes(&self, active_nodes: &[Node]) -> Vec<Node> {
        active_nodes
            .iter()
            .flat_map(|v| self.parents(v))
            // TODO should include active_nodes?
            // .chain(active_nodes.iter())
            .map(|&w| w)
            .unique()
            .collect()
    }
    /// given a list of previously active nodes,
    /// calculates the all candidates of next active nodes
    /// that is all childrens of prev active nodes
    fn f_upconvert_active_nodes<V: VecLike<Prob>>(
        &self,
        prev_layer: &PHMMLayerRaw<V>,
    ) -> Option<Vec<Node>> {
        if let Some(nodes) = &prev_layer.active_nodes {
            let next_active_nodes = self.f_get_next_active_nodes(nodes);
            Some(next_active_nodes)
        } else {
            info!("upconvert=None");
            None
        }
    }
    fn b_upconvert_active_nodes<V: VecLike<Prob>>(
        &self,
        prev_layer: &PHMMLayerRaw<V>,
    ) -> Option<Vec<Node>> {
        if let Some(nodes) = &prev_layer.active_nodes {
            let next_active_nodes = self.b_get_next_active_nodes(nodes);
            Some(next_active_nodes)
        } else {
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
        // get Vec<(index: usize, p: Prob)>
        let mut scores: Vec<(usize, f64)> = match &layer.active_nodes {
            Some(nodes) => nodes
                .iter()
                .map(|v| {
                    let p = layer.pM.get(v.0) + layer.pI.get(v.0) + layer.pD.get(v.0);
                    (v.0, p.to_log_value())
                })
                .collect(),
            None => layer
                .to_kmer_prob()
                .into_iter()
                .map(|p| p.to_log_value())
                .enumerate()
                .collect(),
        };
        // sort descending(from largest to smallest) order
        // TODO sort the vector implement Ord for Prob
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
    fn to_sparse_layer(
        &self,
        param: &PHMMParams,
        dense_layer: &PHMMLayerRaw<DenseVec<Prob>>,
    ) -> PHMMLayerRaw<SparseVec<Prob>> {
        let active_nodes = self
            .refine_active_nodes(dense_layer, param.n_max_active_nodes)
            .unwrap();
        let mut layer = PHMMLayerRaw::<SparseVec<Prob>>::new(self.n_nodes());
        for v in &active_nodes {
            layer.pM.set(v.0, dense_layer.pM.get(v.0));
            layer.pI.set(v.0, dense_layer.pI.get(v.0));
            layer.pD.set(v.0, dense_layer.pD.get(v.0));
        }
        layer.pMB = dense_layer.pMB;
        layer.pIB = dense_layer.pIB;
        layer.pE = dense_layer.pE;
        layer.active_nodes = Some(active_nodes);
        layer
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
        if param.only_active_nodes && index >= param.n_ignore_active_nodes_first {
            next_layer.active_nodes = self.f_upconvert_active_nodes(&prev_layer);
            info!(
                "f_step({}) upconverted {:?}",
                index, next_layer.active_nodes
            );
        } else {
            info!("f_step({}) skipped upconvert", index);
        }
        self.fmi_from_fmid(param, &prev_layer, &mut next_layer, emission);
        self.fb_from_fb(param, &prev_layer, &mut next_layer);
        self.fd_from_fmi(param, &mut next_layer);
        self.fe_from_fmid(param, &mut next_layer);
        // reduce active nodes
        if param.only_active_nodes && index >= param.n_ignore_active_nodes_first {
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
        let mut layers = Vec::new();
        let layer = self.f_init::<DenseVec<Prob>>(&param);
        layers.push(layer);
        for (i, &emission) in emissions.iter().enumerate() {
            let layer = self.f_step(&param, &layers[i], emission, i);
            layers.push(layer);
        }
        layers
    }
    fn forward_sparse(
        &self,
        param: &PHMMParams,
        emissions: &[u8],
    ) -> (
        Vec<PHMMLayerRaw<DenseVec<Prob>>>,
        Vec<PHMMLayerRaw<SparseVec<Prob>>>,
    ) {
        // dense calculation
        let mut layers_first: Vec<PHMMLayerRaw<DenseVec<Prob>>> = Vec::new();
        let layer = self.f_init::<DenseVec<Prob>>(&param);
        layers_first.push(layer);
        let mut iter = emissions.iter().enumerate();
        for (i, &emission) in iter.by_ref().take(param.n_ignore_active_nodes_first) {
            warn!("dense i={}", i);
            let layer = self.f_step(&param, &layers_first[i], emission, i);
            layers_first.push(layer);
        }

        // sparse calculation
        let mut layers: Vec<PHMMLayerRaw<SparseVec<Prob>>> = Vec::new();
        layers.push(self.to_sparse_layer(param, &layers_first[layers_first.len() - 1]));
        for (i, &emission) in iter {
            warn!("sparse i={}", i);
            let layer = self.f_step(
                &param,
                &layers[i - param.n_ignore_active_nodes_first],
                emission,
                i,
            );
            warn!("pE({})={}", i, layer.pE);
            layers.push(layer);
        }
        // remove first element
        layers.remove(0);

        (layers_first, layers)
    }
    /// calculate P(X) i.e. the probability that the PHMM produces the emissions
    /// TODO should be equal to backward_prob
    fn forward_prob(&self, param: &PHMMParams, emissions: &[u8]) -> Prob {
        let layers = self.forward(param, emissions);
        let last_layer = layers.last().unwrap();
        last_layer.pE
    }

    // backward
    /// B[t-1].pD from B[t]
    fn bd_from_bmi<V: VecLike<Prob>>(
        &self,
        param: &PHMMParams,
        layer_t: &PHMMLayerRaw<V>,
        layer_t_1: &mut PHMMLayerRaw<V>,
        emission: u8,
    ) {
        let mut bD0 = V::new(self.n_nodes(), Prob::from_prob(0.0));
        let mut bD1;
        let iterator: Box<dyn std::iter::Iterator<Item = Node>> =
            if let Some(active_nodes) = &layer_t_1.active_nodes {
                Box::new(active_nodes.iter().map(|&v| v))
            } else {
                Box::new(iter_nodes(self.n_nodes()))
            };
        for v in iterator {
            let to_m: Prob = self
                .childs(&v)
                .iter()
                .map(|w| {
                    let emission_prob_M: Prob = if self.emission(w) == emission {
                        param.p_match
                    } else {
                        param.p_mismatch
                    };
                    self.trans_prob(&v, w) * param.p_DM * layer_t.pM.get(w.0) * emission_prob_M
                })
                .sum();
            let emission_prob_I: Prob = param.p_random;
            let to_i = param.p_DI * emission_prob_I * layer_t.pI.get(v.0);
            let bD0v = to_m + to_i;
            bD0.set(v.0, bD0v);
            layer_t_1.pD.set(v.0, bD0v);
        }
        // >0
        for _x in 1..param.n_max_gaps {
            bD1 = V::new(self.n_nodes(), Prob::from_prob(0.0));
            let iterator: Box<dyn std::iter::Iterator<Item = Node>> =
                if let Some(active_nodes) = &layer_t_1.active_nodes {
                    Box::new(active_nodes.iter().map(|&v| v))
                } else {
                    Box::new(iter_nodes(self.n_nodes()))
                };
            for v in iterator {
                let bD1v = self
                    .childs(&v)
                    .iter()
                    .map(|w| self.trans_prob(&v, w) * param.p_DD * bD0.get(w.0))
                    .sum();
                bD1.set(v.0, bD1v);
                layer_t_1.pD.set(v.0, layer_t_1.pD.get(v.0) + bD1v);
            }
            bD0 = bD1;
        }
    }
    /// calc B[t] -> B[t-1]
    fn bmi_from_bmid<V: VecLike<Prob>>(
        &self,
        param: &PHMMParams,
        layer_t: &PHMMLayerRaw<V>,
        layer_t_1: &mut PHMMLayerRaw<V>,
        emission: u8,
    ) {
        let iterator: Box<dyn std::iter::Iterator<Item = Node>> =
            if let Some(active_nodes) = &layer_t_1.active_nodes {
                Box::new(active_nodes.iter().map(|&v| v))
            } else {
                Box::new(iter_nodes(self.n_nodes()))
            };
        for v in iterator {
            // M state
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
                        * (param.p_MM * emission_prob_M * layer_t.pM.get(w.0)
                            + param.p_MD * layer_t_1.pD.get(w.0))
                })
                .sum();
            let emission_prob_I: Prob = param.p_random;
            let to_i = param.p_MI * emission_prob_I * layer_t.pI.get(v.0);
            layer_t_1.pM.set(v.0, to_md + to_i);

            // I state
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
                        * (param.p_IM * emission_prob_M * layer_t.pM.get(w.0)
                            + param.p_ID * layer_t_1.pD.get(w.0))
                })
                .sum();
            let emission_prob_I: Prob = param.p_random;
            let to_i = param.p_II * emission_prob_I * layer_t.pI.get(v.0);
            layer_t_1.pI.set(v.0, to_md + to_i);
        }
    }
    fn bb_from_bmd<V: VecLike<Prob>>(
        &self,
        param: &PHMMParams,
        layer_t: &PHMMLayerRaw<V>,
        layer_t_1: &mut PHMMLayerRaw<V>,
        emission: u8,
    ) {
        // pMB
        let iterator: Box<dyn std::iter::Iterator<Item = Node>> =
            if let Some(active_nodes) = &layer_t_1.active_nodes {
                Box::new(active_nodes.iter().map(|&v| v))
            } else {
                Box::new(iter_nodes(self.n_nodes()))
            };
        let to_md: Prob = iterator
            .map(|v| {
                let emission_prob_M: Prob = if self.emission(&v) == emission {
                    param.p_match
                } else {
                    param.p_mismatch
                };
                self.init_prob(&v)
                    * (param.p_MM * emission_prob_M * layer_t.pM.get(v.0)
                        + param.p_MD * layer_t_1.pD.get(v.0))
            })
            .sum();
        let emission_prob_I: Prob = param.p_random;
        let to_i = param.p_MI * emission_prob_I * layer_t.pIB;
        layer_t_1.pMB = to_md + to_i;

        // pIB
        let iterator: Box<dyn std::iter::Iterator<Item = Node>> =
            if let Some(active_nodes) = &layer_t_1.active_nodes {
                Box::new(active_nodes.iter().map(|&v| v))
            } else {
                Box::new(iter_nodes(self.n_nodes()))
            };
        let to_md: Prob = iterator
            .map(|v| {
                let emission_prob_M: Prob = if self.emission(&v) == emission {
                    param.p_match
                } else {
                    param.p_mismatch
                };
                self.init_prob(&v)
                    * (param.p_IM * emission_prob_M * layer_t.pM.get(v.0)
                        + param.p_ID * layer_t_1.pD.get(v.0))
            })
            .sum();
        let emission_prob_I: Prob = param.p_random;
        let to_i = param.p_II * emission_prob_I * layer_t.pIB;
        layer_t_1.pIB = to_md + to_i;
    }
    fn b_init<V: VecLike<Prob>>(&self, param: &PHMMParams) -> PHMMLayerRaw<V> {
        PHMMLayerRaw::b_init(self.n_nodes(), param.p_end)
    }
    fn b_step<V: VecLike<Prob>>(
        &self,
        param: &PHMMParams,
        prev_layer: &PHMMLayerRaw<V>,
        emission: u8,
        index: usize,
    ) -> PHMMLayerRaw<V> {
        // emission: x[i+1]
        let mut next_layer = PHMMLayerRaw::new(self.n_nodes());
        if param.only_active_nodes && index >= param.n_ignore_active_nodes_first {
            next_layer.active_nodes = self.b_upconvert_active_nodes(&prev_layer);
            info!(
                "b_step({}) upconverted {:?}",
                index, next_layer.active_nodes
            );
        } else {
            info!("b_step({}) skipped upconvert", index);
        }
        self.bd_from_bmi(param, &prev_layer, &mut next_layer, emission);
        self.bmi_from_bmid(param, &prev_layer, &mut next_layer, emission);
        self.bb_from_bmd(param, &prev_layer, &mut next_layer, emission);
        if param.only_active_nodes && index >= param.n_ignore_active_nodes_first {
            next_layer.active_nodes =
                self.refine_active_nodes(&next_layer, param.n_max_active_nodes);
            info!("b_step({}) refined {:?}", index, next_layer.active_nodes);
            // FIXME
            // next_layer.active_nodes = None
        }
        // check the cache miss rate
        next_layer
    }
    /// Backward probability
    /// layers[i] = P(x[i..] | start with a state)
    fn backward(&self, param: &PHMMParams, emissions: &[u8]) -> Vec<PHMMLayer> {
        debug!("start backward!");
        let mut layers = Vec::new();
        let layer = self.b_init::<DenseVec<Prob>>(&param);
        trace!("l0:{}", layer.pM.len());
        trace!("l0:\n{}", layer);
        layers.push(layer);
        for (i, &emission) in emissions.iter().rev().enumerate() {
            info!("l{} emission={}", i, emission as char);
            let layer = self.b_step(&param, &layers[i], emission, i);
            trace!("l{}:\n{}", i, layer);
            layers.push(layer);
        }
        layers.reverse();
        layers
    }
    fn backward_sparse(
        &self,
        param: &PHMMParams,
        emissions: &[u8],
    ) -> (
        Vec<PHMMLayerRaw<SparseVec<Prob>>>,
        Vec<PHMMLayerRaw<DenseVec<Prob>>>,
    ) {
        // dense
        let mut layers_last: Vec<PHMMLayerRaw<DenseVec<Prob>>> = Vec::new();
        let layer = self.b_init::<DenseVec<Prob>>(&param);
        layers_last.push(layer);
        let mut iter = emissions.iter().rev().enumerate();
        for (i, &emission) in iter.by_ref().take(param.n_ignore_active_nodes_first) {
            warn!("dense i={}", i);
            let layer = self.b_step(&param, &layers_last[i], emission, i);
            layers_last.push(layer);
        }

        // sparse
        let mut layers: Vec<PHMMLayerRaw<SparseVec<Prob>>> = Vec::new();
        layers.push(self.to_sparse_layer(param, &layers_last[layers_last.len() - 1]));
        for (i, &emission) in iter {
            warn!("sparse i={}", i);
            let layer = self.b_step(
                &param,
                &layers[i - param.n_ignore_active_nodes_first],
                emission,
                i,
            );
            layers.push(layer);
        }
        // remove first element
        layers.remove(0);

        // reverse the vector
        layers.reverse();
        layers_last.reverse();

        (layers, layers_last)
    }
    /// calculate P(X) i.e. the probability that the PHMM produces the emissions
    /// TODO should be equal to forward_prob
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
}

#[cfg(test)]
mod tests {
    use super::*;
}
