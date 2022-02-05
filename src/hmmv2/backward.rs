//!
//! Backward algorithm definitions
//!

use super::common::{PHMMEdge, PHMMModel, PHMMNode};
use super::table::{PHMMResult, PHMMTable};
use crate::prob::Prob;
use crate::vector::{NodeVec, Storage};

///
/// Backward Algorithm
///
impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    /// Run Backward algorithm to the emissions
    ///
    pub fn backward<S>(&self, emissions: &[u8]) -> PHMMResult<S>
    where
        S: Storage<Item = Prob>,
    {
        let r0 = PHMMResult {
            init_table: self.b_init(),
            tables: Vec::new(),
        };
        // feed the emissions backward
        let mut r = emissions
            .iter()
            .rev()
            .enumerate()
            .fold(r0, |mut r, (i, &emission)| {
                let table = if i == 0 {
                    self.b_step(i, emission, &r.init_table)
                } else {
                    self.b_step(i, emission, r.tables.last().unwrap())
                };
                r.tables.push(table);
                r
            });
        // reverse the vector, to order the tables along with emissions
        // i.e. tables[i] corresponds to the emissions[i]
        r.tables.reverse();
        r
    }
    ///
    /// Create init_table in PHMMResult for Backward algorithm
    ///
    /// ```text
    /// b*_n[k]
    /// = P(go to end | starts from state *_k)
    /// = p_e (if *_k is reachable to end)
    ///     or
    ///   0   (otherwise)
    /// ```
    ///
    fn b_init<S>(&self) -> PHMMTable<S>
    where
        S: Storage<Item = Prob>,
    {
        let p_end = self.param.p_end;
        PHMMTable::new(
            self.n_nodes(),
            // m,i,d is p_end
            p_end,
            p_end,
            p_end,
            // begin,end is 0
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
        )
    }
    ///
    /// Calculate the table from the previous table
    /// for Backward algorithm
    ///
    fn b_step<S>(&self, i: usize, emission: u8, prev_table: &PHMMTable<S>) -> PHMMTable<S>
    where
        S: Storage<Item = Prob>,
    {
        let mut table = PHMMTable::new(
            self.n_nodes(),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
        );
        // bd (silent states) should be first
        self.bd(&mut table, prev_table, emission);
        self.be(&mut table, prev_table, emission);
        // normal state is next
        self.bm(&mut table, prev_table, emission);
        self.bi(&mut table, prev_table, emission);
        self.bib(&mut table, prev_table, emission);
        self.bmb(&mut table, prev_table, emission);
        table
    }
}

// functions to calculate each step
impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    /// Fill the backward probs of `Del` states
    ///
    /// For `i=0,...,n-1`
    ///
    /// ```text
    /// bd_i[v]
    /// = P(emits x[i:] | starts from state d_v)
    /// =   (to_m_childs) \sum_{w: childs} t_vw p_dm e(x[i]) bm_i+1[w]
    ///   + (to_d_childs) \sum_{w: childs} t_vw p_dd         bd_i[w]
    ///   + (to_i_self)                         p_di e(x[i]) bi_i+1[v]
    /// ```
    /// (Here `x[i:] = x[i],...,x[n-1]`)
    ///
    /// This calculation has a recursive definition of `bd_i`, so purging them
    /// by allowing only `param.n_max_gaps` times continuous deletions.
    ///
    /// ```text
    /// bd_i(0)[v]
    /// =   (to_m_childs) \sum_{w: childs} t_vw p_dm e(x[i]) bm_i+1[w]
    ///   + (to_i_self)                         p_di e(x[i]) bi_i+1[v]
    ///
    /// t = 1,2,3,...
    /// bd_i(t)[v]
    /// =   (to_d_childs) \sum_{w: childs} t_vw p_dd         bd_i(t-1)[w]
    ///
    /// bd_i[v] = \sum_t bd_i(t)[v]
    /// ```
    ///
    /// ## Dependency
    ///
    /// It depends on `bm_i+1, bi_i+1`. This can be calculated first in `i` tables.
    ///
    fn bd<S>(&self, t0: &mut PHMMTable<S>, t1: &PHMMTable<S>, emission: u8)
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        let mut bdt0 = self.bd0(t0, emission);
        t0.d += &bdt0;
        for t in 0..param.n_max_gaps {
            bdt0 = self.bdt(&bdt0);
            t0.d += &bdt0;
        }
    }
    ///
    /// Calculate `bd_i(t=0)[v]`
    ///
    /// ```text
    /// bd_i(0)[v]
    /// =   (to_m_childs) \sum_{w: childs} t_vw p_dm e(x[i]) bm_i+1[w]
    ///   + (to_i_self)                         p_di e(x[i]) bi_i+1[v]
    /// ```
    ///
    /// (For the complete definitions, see the reference doc of `self.bd()`)
    fn bd0<S>(&self, t0: &PHMMTable<S>, emission: u8) -> NodeVec<S>
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        let mut bd0 = NodeVec::new(self.n_nodes(), Prob::from_prob(0.0));
        for (k, kw) in self.nodes() {
            // (1) to match
            let p_to_match: Prob = self
                .childs(k)
                .map(|(_, l, ew)| {
                    // k -> l
                    let p_trans = ew.trans_prob();
                    let lw = self.graph.node_weight(l).unwrap();
                    // emission prob on l
                    let p_emit = if lw.emission() == emission {
                        param.p_match
                    } else {
                        param.p_mismatch
                    };
                    p_trans * param.p_DM * p_emit * t0.m[l]
                })
                .sum();

            // (2) to ins
            let p_to_ins = param.p_DI * param.p_random * t0.i[k];
            bd0[k] = p_to_match + p_to_ins;
        }
        bd0
    }
    ///
    /// Calculate `bd_i(t)[v]` for `t > 0`
    ///
    /// ```text
    /// bd_i(t)[v]
    /// =   (to_d_childs) \sum_{w: childs} t_vw p_dd         bd_i(t-1)[w]
    /// ```
    ///
    /// (For the complete definitions, see the reference doc of `self.bd()`)
    fn bdt<S>(&self, bdt1: &NodeVec<S>) -> NodeVec<S>
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        let mut bdt0 = NodeVec::new(self.n_nodes(), Prob::from_prob(0.0));
        for (k, kw) in self.nodes() {
            bdt0[k] = self
                .childs(k)
                .map(|(_, l, ew)| {
                    // k -> l
                    let p_trans = ew.trans_prob();
                    p_trans * param.p_DD * bdt1[l]
                })
                .sum();
        }
        bdt0
    }
    /// Fill the backward probs of `Match` states
    ///
    /// For `i=0,...,n-1`
    ///
    /// ```text
    /// bm_i[v]
    /// = P(emits x[i:] | starts from state m_v)
    /// =   (to_m_childs) \sum_{w: childs} t_vw p_mm e(x[i]) bm_i+1[w]
    ///   + (to_d_childs) \sum_{w: childs} t_vw p_md         bd_i[w]
    ///   + (to_i_self)                         p_mi e(x[i]) bi_i+1[v]
    /// ```
    ///
    /// (Here `x[i:] = x[i],...,x[n-1]`)
    ///
    /// ## Dependency
    ///
    /// It depends on `bm_i+1, bi_i+1, bd_i`.
    ///
    fn bm<S>(&self, t0: &mut PHMMTable<S>, t1: &PHMMTable<S>, emission: u8)
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        for (k, kw) in self.nodes() {
            // (1) to match and del
            let p_to_match_del: Prob = self
                .childs(k)
                .map(|(_, l, ew)| {
                    // k -> l
                    let p_trans = ew.trans_prob();
                    let lw = self.graph.node_weight(l).unwrap();
                    // emission prob on l
                    let p_emit = if lw.emission() == emission {
                        param.p_match
                    } else {
                        param.p_mismatch
                    };
                    p_trans * ((param.p_MM * p_emit * t1.m[l]) + (param.p_MD * t0.d[l]))
                })
                .sum();

            // (2) to ins
            let p_to_ins = param.p_MI * param.p_random * t1.i[k];

            // sum
            t0.m[k] = p_to_match_del + p_to_ins;
        }
    }
    /// Fill the backward probs of `Ins` states
    ///
    /// For `i=0,...,n-1`
    ///
    /// ```text
    /// bi_i[v]
    /// = P(emits x[i:] | starts from state i_v)
    /// =   (to_m_childs) \sum_{w: childs} t_vw p_im e(x[i]) bm_i+1[w]
    ///   + (to_d_childs) \sum_{w: childs} t_vw p_id         bd_i[w]
    ///   + (to_i_self)                         p_ii e(x[i]) bi_i+1[v]
    /// ```
    /// (Here `x[i:] = x[i],...,x[n-1]`)
    ///
    /// ## Dependency
    ///
    /// It depends on `bm_i+1, bi_i+1, bd_i`.
    ///
    fn bi<S>(&self, t0: &mut PHMMTable<S>, t1: &PHMMTable<S>, emission: u8)
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        for (k, kw) in self.nodes() {
            // (1) to match and del
            let p_to_match_del: Prob = self
                .childs(k)
                .map(|(_, l, ew)| {
                    // k -> l
                    let p_trans = ew.trans_prob();
                    let lw = self.graph.node_weight(l).unwrap();
                    // emission prob on l
                    let p_emit = if lw.emission() == emission {
                        param.p_match
                    } else {
                        param.p_mismatch
                    };
                    p_trans * ((param.p_IM * p_emit * t1.m[l]) + (param.p_ID * t0.d[l]))
                })
                .sum();

            // (2) to ins
            let p_to_ins = param.p_II * param.p_random * t1.i[k];

            // sum
            t0.i[k] = p_to_match_del + p_to_ins;
        }
    }
    /// fill the backward prob of `MatchBegin` state
    ///
    /// ```text
    /// bm_i[b]
    /// = P(emits x[i:] | starts from state m_b (MatchBegin))
    /// =   (to_all_m) \sum_{w} t_bw p_mm e(x[i]) bm_i+1[w]
    ///   + (to_all_d) \sum_{w} t_bw p_md         bd_i[w]
    ///   + (to_self_i)              p_mi e(x[i]) bi_i+1[b]
    /// ```
    ///
    /// Here `t_bw` is a init probability from Begin to node `w`.
    ///
    /// ## Dependency
    /// bm_i+1, bi_i+1, bd_i
    ///
    fn bmb<S>(&self, t0: &mut PHMMTable<S>, t1: &PHMMTable<S>, emission: u8)
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        for (k, kw) in self.nodes() {
            // (1) to match and del of all nodes
            let p_to_match_del: Prob = self
                .nodes()
                .map(|(l, lw)| {
                    // k=Begin -> l
                    let p_trans = lw.init_prob();
                    // emission prob on l
                    let p_emit = if lw.emission() == emission {
                        param.p_match
                    } else {
                        param.p_mismatch
                    };
                    p_trans * ((param.p_MM * p_emit * t1.m[l]) + (param.p_MD * t0.d[l]))
                })
                .sum();

            // (2) to ins of self node
            let p_to_ins = param.p_MI * param.p_random * t1.ib;

            // sum
            t0.mb = p_to_match_del + p_to_ins;
        }
    }
    /// fill the backward prob of `InsBegin` state
    ///
    /// ```text
    /// bi_i[b]
    /// = P(emits x[i:] | starts from state i_b (InsBegin))
    /// =   (to_all_m) \sum_{w} t_bw p_im e(x[i]) bm_i+1[w]
    ///   + (to_all_d) \sum_{w} t_bw p_id         bd_i[w]
    ///   + (to_self_i)              p_ii e(x[i]) bi_i+1[b]
    /// ```
    ///
    /// Here `t_bw` is a init probability from Begin to node `w`.
    ///
    /// ## Dependency
    /// bm_i+1, bi_i+1, bd_i
    ///
    fn bib<S>(&self, t0: &mut PHMMTable<S>, t1: &PHMMTable<S>, emission: u8)
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        for (k, kw) in self.nodes() {
            // (1) to match and del of all nodes
            let p_to_match_del: Prob = self
                .nodes()
                .map(|(l, lw)| {
                    // k=Begin -> l
                    let p_trans = lw.init_prob();
                    // emission prob on l
                    let p_emit = if lw.emission() == emission {
                        param.p_match
                    } else {
                        param.p_mismatch
                    };
                    p_trans * ((param.p_IM * p_emit * t1.m[l]) + (param.p_ID * t0.d[l]))
                })
                .sum();

            // (2) to ins of self node
            let p_to_ins = param.p_II * param.p_random * t1.ib;

            // sum
            t0.ib = p_to_match_del + p_to_ins;
        }
    }
    /// Fill the backward prob of `End` state
    ///
    /// `be_i[v] = 0` because `End` has no emission and no childs
    ///
    /// ## Dependency
    /// None
    ///
    fn be<S>(&self, t0: &mut PHMMTable<S>, t1: &PHMMTable<S>, _emission: u8)
    where
        S: Storage<Item = Prob>,
    {
        t0.e = Prob::from_prob(0.0);
    }
}

#[cfg(test)]
mod tests {
    use super::super::seqgraph::create_linear_seq_graph;
    use super::*;

    #[test]
    fn create_linear_seq_graph_test() {
        let g = create_linear_seq_graph(b"ATCGGCTAGC");
        let phmm = g.to_phmm();
        println!("{}", phmm);
    }
}
