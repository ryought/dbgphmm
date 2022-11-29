//!
//! Backward algorithm definitions
//!

use super::common::{PHMMEdge, PHMMModel, PHMMNode};
use super::hint::Hint;
use super::result::{PHMMResult, PHMMResultLike, PHMMResultSparse};
use super::table::PHMMTable;
use crate::graph::active_nodes::ActiveNodes;
use crate::prob::Prob;
use crate::vector::{NodeVec, Storage};

///
/// Backward Algorithm
///
impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    /// Run Backward algorithm to the emissions
    ///
    /// `bt_i[k]` = P(emits `x[i:] = x[i], ..., x[n-1]` | starts from state `t_k`)
    ///
    /// * `t` is a type of state, either Match, Ins, Del
    /// * `k` is a node index
    ///
    pub fn backward(&self, emissions: &[u8]) -> PHMMResult {
        let r0 = PHMMResult {
            init_table: self.b_init(),
            tables: Vec::new(),
            is_forward: false,
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
    pub fn backward_with_hint(&self, emissions: &[u8], hint: &Hint) -> PHMMResult {
        let r0 = PHMMResult {
            init_table: self.b_init(),
            tables: Vec::new(),
            is_forward: false,
        };
        let n = emissions.len();
        // feed the emissions backward
        let mut r = emissions
            .iter()
            .enumerate()
            .rev()
            .fold(r0, |mut r, (i, &emission)| {
                let table = if i == n - 1 {
                    self.b_step_with_active_nodes(i, emission, &r.init_table, hint.active_nodes(i))
                } else {
                    self.b_step_with_active_nodes(
                        i,
                        emission,
                        r.tables.last().unwrap(),
                        hint.active_nodes(i),
                    )
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
    /// Run Backward algorithm to the emissions, with sparse calculation
    ///
    pub fn backward_sparse(&self, emissions: &[u8]) -> PHMMResultSparse {
        let r0 = PHMMResultSparse {
            init_table: self.b_init(),
            tables_warmup: Vec::new(),
            tables_sparse: Vec::new(),
            is_forward: false,
        };
        let param = &self.param;
        let mut r = emissions
            .iter()
            .rev()
            .enumerate()
            .fold(r0, |mut r, (i, &emission)| {
                if i < param.n_warmup {
                    // dense_table -> dense_table
                    let table = if i == 0 {
                        self.b_step(i, emission, &r.init_table)
                    } else {
                        self.b_step(i, emission, r.tables_warmup.last().unwrap())
                    };
                    r.tables_warmup.push(table);
                } else if i == param.n_warmup {
                    // dense_table -> sparse_table
                    let table_prev = r
                        .tables_warmup
                        .last()
                        .unwrap()
                        .to_sparse_active_nodes(param.n_active_nodes);
                    let mut table = self.b_step(i, emission, &table_prev);
                    table.refresh_active_nodes(param.n_active_nodes);
                    r.tables_sparse.push(table);
                } else {
                    // sparse_table -> sparse_table
                    let mut table = self.b_step(i, emission, r.tables_sparse.last().unwrap());
                    table.refresh_active_nodes(param.n_active_nodes);
                    r.tables_sparse.push(table);
                };
                r
            });
        // reverse the vector, to order the tables along with emissions
        // i.e. tables[i] corresponds to the emissions[i]
        r.tables_warmup.reverse();
        r.tables_sparse.reverse();
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
    fn b_step<S>(&self, _i: usize, emission: u8, prev_table: &PHMMTable<S>) -> PHMMTable<S>
    where
        S: Storage<Item = Prob>,
    {
        // active_nodes are not used in bd and be
        let mut table = PHMMTable::new(
            self.n_nodes(),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
        );
        let active_nodes = &prev_table.active_nodes;

        // bd (silent states) should be first
        self.bd(&mut table, prev_table, emission, active_nodes);
        self.be(&mut table, prev_table, emission, active_nodes);

        // candidates of active nodes of next step
        let active_nodes = if !S::is_dense() {
            let parents_and_us = prev_table.active_nodes.to_parents_and_us(self);
            let active_in_d = table.active_nodes_from_prob(self.param.n_active_nodes);
            parents_and_us.merge(&active_in_d)
        } else {
            ActiveNodes::All
        };

        // normal state is next
        self.bm(&mut table, prev_table, emission, &active_nodes);
        self.bi(&mut table, prev_table, emission, &active_nodes);
        self.bib(&mut table, prev_table, emission, &active_nodes);
        self.bmb(&mut table, prev_table, emission, &active_nodes);
        table
    }
    ///
    /// Calculate the table from the previous table
    /// for Backward algorithm
    ///
    fn b_step_with_active_nodes<S>(
        &self,
        _i: usize,
        emission: u8,
        prev_table: &PHMMTable<S>,
        active_nodes: &ActiveNodes,
    ) -> PHMMTable<S>
    where
        S: Storage<Item = Prob>,
    {
        // active_nodes are not used in bd and be
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
        self.bd_with_active_nodes(&mut table, prev_table, emission, &active_nodes);
        self.be(&mut table, prev_table, emission, &active_nodes);

        // normal state is next
        self.bm(&mut table, prev_table, emission, &active_nodes);
        self.bi(&mut table, prev_table, emission, &active_nodes);
        self.bib(&mut table, prev_table, emission, &active_nodes);
        self.bmb(&mut table, prev_table, emission, &active_nodes);
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
    fn bd<S>(
        &self,
        t0: &mut PHMMTable<S>,
        t1: &PHMMTable<S>,
        emission: u8,
        active_nodes: &ActiveNodes,
    ) where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;

        // bd0.d[k] depends on child and itself.
        let mut active_nodes_t = active_nodes.to_parents_and_us(self);
        let mut bdt0 = self.bd0(t1, emission, &active_nodes_t);
        t0.d += &bdt0.d;

        for _t in 0..param.n_max_gaps {
            // bdt0.d[k] only depends on k's child l in bdt1.d
            if active_nodes_t.is_only() {
                active_nodes_t = bdt0.active_nodes_from_prob(param.n_active_nodes);
                active_nodes_t = active_nodes_t.to_parents(self);
            };
            bdt0 = self.bdt(&bdt0, &active_nodes_t);
            t0.d += &bdt0.d;
        }
    }
    ///
    /// `bd` with fixed active nodes (not growing)
    ///
    fn bd_with_active_nodes<S>(
        &self,
        t0: &mut PHMMTable<S>,
        t1: &PHMMTable<S>,
        emission: u8,
        active_nodes: &ActiveNodes,
    ) where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;

        // bd0.d[k] depends on child and itself.
        let mut bdt0 = self.bd0(t1, emission, &active_nodes);
        t0.d += &bdt0.d;

        for _t in 0..param.n_max_gaps {
            // bdt0.d[k] only depends on k's child l in bdt1.d
            bdt0 = self.bdt(&bdt0, &active_nodes);
            t0.d += &bdt0.d;
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
    fn bd0<S>(&self, t0: &PHMMTable<S>, emission: u8, active_nodes: &ActiveNodes) -> PHMMTable<S>
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;

        let mut bd0 = PHMMTable::zero(self.n_nodes());
        for (k, _) in self.active_nodes(active_nodes) {
            // (1) to match
            let p_to_match: Prob = self
                .childs(k)
                .map(|(_, l, ew)| {
                    // k -> l
                    let p_trans = ew.trans_prob();
                    // emission prob on l
                    let p_emit = self.p_match_emit(l, emission);
                    p_trans * param.p_DM * p_emit * t0.m[l]
                })
                .sum();

            // (2) to ins
            let p_to_ins = param.p_DI * self.p_ins_emit() * t0.i[k];
            bd0.d[k] = p_to_match + p_to_ins;
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
    fn bdt<S>(&self, bdt1: &PHMMTable<S>, active_nodes: &ActiveNodes) -> PHMMTable<S>
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;

        let mut bdt0 = PHMMTable::zero(self.n_nodes());

        for (k, _) in self.active_nodes(active_nodes) {
            bdt0.d[k] = self
                .childs(k)
                .map(|(_, l, ew)| {
                    // k -> l
                    let p_trans = ew.trans_prob();
                    p_trans * param.p_DD * bdt1.d[l]
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
    fn bm<S>(
        &self,
        t0: &mut PHMMTable<S>,
        t1: &PHMMTable<S>,
        emission: u8,
        active_nodes: &ActiveNodes,
    ) where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        for (k, _) in self.active_nodes(active_nodes) {
            // (1) to match and del
            let p_to_match_del: Prob = self
                .childs(k)
                .map(|(_, l, ew)| {
                    // k -> l
                    let p_trans = ew.trans_prob();
                    // emission prob on l
                    let p_emit = self.p_match_emit(l, emission);
                    p_trans * ((param.p_MM * p_emit * t1.m[l]) + (param.p_MD * t0.d[l]))
                })
                .sum();

            // (2) to ins
            let p_to_ins = param.p_MI * self.p_ins_emit() * t1.i[k];

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
    fn bi<S>(
        &self,
        t0: &mut PHMMTable<S>,
        t1: &PHMMTable<S>,
        emission: u8,
        active_nodes: &ActiveNodes,
    ) where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        for (k, _) in self.active_nodes(active_nodes) {
            // (1) to match and del
            let p_to_match_del: Prob = self
                .childs(k)
                .map(|(_, l, ew)| {
                    // k -> l
                    let p_trans = ew.trans_prob();
                    // emission prob on l
                    let p_emit = self.p_match_emit(l, emission);
                    p_trans * ((param.p_IM * p_emit * t1.m[l]) + (param.p_ID * t0.d[l]))
                })
                .sum();

            // (2) to ins
            let p_to_ins = param.p_II * self.p_ins_emit() * t1.i[k];

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
    fn bmb<S>(
        &self,
        t0: &mut PHMMTable<S>,
        t1: &PHMMTable<S>,
        emission: u8,
        active_nodes: &ActiveNodes,
    ) where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;

        // (1) to match and del of all nodes
        let p_to_match_del: Prob = self
            .active_nodes(active_nodes)
            .map(|(l, lw)| {
                // k=Begin -> l
                let p_trans = lw.init_prob();
                // emission prob on l
                let p_emit = self.p_match_emit(l, emission);
                p_trans * ((param.p_MM * p_emit * t1.m[l]) + (param.p_MD * t0.d[l]))
            })
            .sum();

        // (2) to ins of self node
        let p_to_ins = param.p_MI * self.p_ins_emit() * t1.ib;

        // sum
        t0.mb = p_to_match_del + p_to_ins;
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
    fn bib<S>(
        &self,
        t0: &mut PHMMTable<S>,
        t1: &PHMMTable<S>,
        emission: u8,
        active_nodes: &ActiveNodes,
    ) where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;

        // (1) to match and del of all nodes
        let p_to_match_del: Prob = self
            .active_nodes(&active_nodes)
            .map(|(l, lw)| {
                // k=Begin -> l
                let p_trans = lw.init_prob();
                // emission prob on l
                let p_emit = self.p_match_emit(l, emission);
                p_trans * ((param.p_IM * p_emit * t1.m[l]) + (param.p_ID * t0.d[l]))
            })
            .sum();

        // (2) to ins of self node
        let p_to_ins = param.p_II * self.p_ins_emit() * t1.ib;

        // sum
        t0.ib = p_to_match_del + p_to_ins;
    }
    /// Fill the backward prob of `End` state
    ///
    /// `be_i[v] = 0` because `End` has no emission and no childs
    ///
    /// ## Dependency
    /// None
    ///
    fn be<S>(
        &self,
        t0: &mut PHMMTable<S>,
        _t1: &PHMMTable<S>,
        _emission: u8,
        _active_nodes: &ActiveNodes,
    ) where
        S: Storage<Item = Prob>,
    {
        t0.e = Prob::from_prob(0.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::ni;
    use crate::hmmv2::mocks::mock_linear_phmm;
    use crate::hmmv2::params::PHMMParams;
    use crate::prob::lp;
    use crate::vector::DenseStorage;
    #[test]
    fn hmm_backward_mock_linear_zero_error() {
        let params = PHMMParams::zero_error();
        println!("{}", params);
        let phmm = mock_linear_phmm(params);
        let r = phmm.backward(b"CGATC");
        for table in r.tables.iter() {
            println!("{}", table);
        }
        println!("{}", r.init_table);
        // total probability
        assert_abs_diff_eq!(r.tables[0].mb, lp(-13.8155605), epsilon = 0.00001);
        // position-wise
        assert_abs_diff_eq!(r.tables[4].m[ni(6)], lp(-11.5129354), epsilon = 0.00001);
        assert_abs_diff_eq!(r.tables[4].m[ni(2)], lp(-11.5129354), epsilon = 0.00001);
        assert_abs_diff_eq!(r.tables[3].m[ni(5)], lp(-11.5129454), epsilon = 0.00001);
        assert_abs_diff_eq!(r.tables[3].m[ni(1)], lp(-11.5129454), epsilon = 0.00001);
        assert_abs_diff_eq!(r.tables[2].m[ni(4)], lp(-11.5129554), epsilon = 0.00001);
        assert_abs_diff_eq!(r.tables[1].m[ni(3)], lp(-11.5129654), epsilon = 0.00001);
        assert_abs_diff_eq!(r.tables[0].m[ni(2)], lp(-11.5129754), epsilon = 0.00001);
        // with allowing no errors, CGATT cannot be emitted.
        // so it should have p=0
        let r2 = phmm.backward(b"CGATT");
        assert_eq!(r2.tables.len(), 5);
        assert!(r2.tables[0].mb.is_zero());
        for table in r2.tables.iter() {
            println!("{}", table);
        }
        println!("{}", r2.init_table);
    }
    #[test]
    fn hmm_backward_mock_linear_high_error() {
        let mut param = PHMMParams::high_error();
        param.n_warmup = 2;
        param.n_active_nodes = 2;
        let phmm = mock_linear_phmm(param);
        // read 1
        let r = phmm.backward(b"CGATC");
        for table in r.tables.iter() {
            println!("{}", table);
        }
        println!("{}", r.init_table);
        assert_abs_diff_eq!(r.tables[0].m[ni(2)], lp(-13.0679200), epsilon = 0.00001);
        assert_abs_diff_eq!(r.tables[0].mb, lp(-15.2115765494), epsilon = 0.00001);
        // read 2
        let r2 = phmm.backward(b"CGATT");
        assert_eq!(r2.tables.len(), 5);
        for table in r2.tables.iter() {
            println!("{}", table);
        }
        println!("{}", r2.init_table);
        assert_abs_diff_eq!(r2.tables[0].mb, lp(-16.7787277), epsilon = 0.00001);
    }
    #[test]
    fn hmm_backward_with_hint_mock_linear_high_error() {
        let phmm = mock_linear_phmm(PHMMParams::high_error());
        // read1
        let read1 = b"CGATC";
        let o = phmm.run(read1);
        let hint = o.to_hint(5);
        println!("{:?}", hint);
        println!("{:?}", hint.len());

        let r1 = phmm.backward(read1);
        let r2 = phmm.backward_with_hint(read1, &hint);
        println!("{}", r1.last_table());
        println!("{}", r2.last_table());
    }
}
