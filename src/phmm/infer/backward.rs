use super::super::PHMM;
use super::table::{PHMMTable, PHMMTables};
use crate::prob::{p, Prob};
use petgraph::graph::NodeIndex;

///
/// Backward algorithm
///
impl PHMM {
    ///
    /// Run Backward algorithm to the emissions
    ///
    /// `bt_i[k]` = P(emits `x[i:] = x[i], ..., x[n-1]` | starts from state `t_k`)
    ///
    /// * `t` is a type of state, either Match, Ins, Del
    /// * `k` is a node index
    ///
    pub fn backward(&self, emissions: &[u8]) -> PHMMTables {
        let r0 = PHMMTables {
            init_table: self.b_init(true),
            tables: Vec::new(),
            is_forward: false,
        };
        let all_nodes = self.to_all_nodes();
        // feed the emissions backward
        let mut r = emissions
            .iter()
            .rev()
            .enumerate()
            .fold(r0, |mut r, (i, &emission)| {
                let table_prev = if i == 0 {
                    &r.init_table
                } else {
                    r.last_table()
                };
                let table = self.b_step(i, emission, table_prev, &all_nodes, true, false);
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
    pub fn backward_sparse(&self, emissions: &[u8]) -> PHMMTables {
        let r0 = PHMMTables {
            init_table: self.b_init(true),
            tables: Vec::new(),
            is_forward: false,
        };
        let param = &self.param;
        let all_nodes = self.to_all_nodes();
        let mut r = emissions
            .iter()
            .rev()
            .enumerate()
            .fold(r0, |mut r, (i, &emission)| {
                if i < param.n_warmup {
                    // dense_table
                    let table_prev = if i == 0 {
                        &r.init_table
                    } else {
                        r.last_table()
                    };
                    let table = self.b_step(i, emission, table_prev, &all_nodes, true, false);
                    r.tables.push(table);
                } else {
                    // sparse_table
                    let table_prev = r.last_table();
                    let active_nodes = table_prev.top_nodes(param.n_active_nodes);
                    let table = self.b_step(i, emission, &table_prev, &active_nodes, false, true);
                    r.tables.push(table);
                };
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
    fn b_init(&self, is_dense: bool) -> PHMMTable {
        let p_end = self.param.p_end;
        PHMMTable::new(
            is_dense,
            self.n_nodes(),
            // m,i,d is p_end
            p_end,
            p_end,
            p_end,
            // begin,end is 0
            p(0.0),
            p(0.0),
            p(0.0),
        )
    }
    ///
    /// Calculate the table from the previous table
    /// for Backward algorithm
    ///
    fn b_step(
        &self,
        i: usize,
        emission: u8,
        prev_table: &PHMMTable,
        nodes: &[NodeIndex],
        is_dense: bool,
        is_adaptive: bool,
    ) -> PHMMTable {
        // active_nodes are not used in bd and be
        let mut table = PHMMTable::new(
            is_dense,
            self.n_nodes(),
            p(0.0),
            p(0.0),
            p(0.0),
            p(0.0),
            p(0.0),
            p(0.0),
        );

        // bd (silent states) should be first
        self.bd(&mut table, prev_table, emission, nodes, is_adaptive);
        self.be(&mut table, prev_table, emission, nodes);

        // candidates of active nodes of next step
        let mut active_nodes = None;
        if is_adaptive {
            let parents_and_us = self.to_parents_and_us(nodes);
            active_nodes = Some(parents_and_us);
            // let active_in_d = table.active_nodes_from_prob(self.param.n_active_nodes);
            // active_nodes = Some(self.merge(&parents_and_us, active_in_d.nodes()));
        };

        // normal state is next
        let active_nodes_ref = if is_adaptive {
            active_nodes.as_ref().unwrap()
        } else {
            nodes
        };
        self.bm(&mut table, prev_table, emission, &active_nodes_ref);
        self.bi(&mut table, prev_table, emission, &active_nodes_ref);
        self.bib(&mut table, prev_table, emission, &active_nodes_ref);
        self.bmb(&mut table, prev_table, emission, &active_nodes_ref);
        table
    }
}

// functions to calculate each step
impl PHMM {
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
    fn bd(
        &self,
        t0: &mut PHMMTable,
        t1: &PHMMTable,
        emission: u8,
        nodes: &[NodeIndex],
        is_adaptive: bool,
    ) {
        let param = &self.param;
        let mut active_nodes = None;
        let is_dense = t0.is_dense();

        // bd0.d[k] depends on child and itself.
        if is_adaptive {
            active_nodes = Some(self.to_parents_and_us(&nodes));
        }
        let mut bdt0 = self.bd0(
            t1,
            emission,
            if is_adaptive {
                active_nodes.as_ref().unwrap()
            } else {
                nodes
            },
            is_dense,
        );
        t0.d += &bdt0.d;

        for _t in 0..param.n_max_gaps {
            // bdt0.d[k] only depends on k's child l in bdt1.d
            if is_adaptive {
                active_nodes = Some(self.to_parents_and_us(active_nodes.as_ref().unwrap()));
            }
            bdt0 = self.bdt(
                &bdt0,
                if is_adaptive {
                    active_nodes.as_ref().unwrap()
                } else {
                    nodes
                },
                is_dense,
            );
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
    fn bd0(&self, t0: &PHMMTable, emission: u8, nodes: &[NodeIndex], is_dense: bool) -> PHMMTable {
        let param = &self.param;

        let mut bd0 = PHMMTable::zero(is_dense, self.n_nodes());
        for &k in nodes {
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
    fn bdt(&self, bdt1: &PHMMTable, nodes: &[NodeIndex], is_dense: bool) -> PHMMTable {
        let param = &self.param;

        let mut bdt0 = PHMMTable::zero(is_dense, self.n_nodes());

        for &k in nodes {
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
    fn bm(&self, t0: &mut PHMMTable, t1: &PHMMTable, emission: u8, nodes: &[NodeIndex]) {
        let param = &self.param;
        for &k in nodes {
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
    fn bi(&self, t0: &mut PHMMTable, t1: &PHMMTable, emission: u8, nodes: &[NodeIndex]) {
        let param = &self.param;
        for &k in nodes {
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
    fn bmb(&self, t0: &mut PHMMTable, t1: &PHMMTable, emission: u8, nodes: &[NodeIndex]) {
        let param = &self.param;

        // (1) to match and del of all nodes
        let p_to_match_del: Prob = nodes
            .iter()
            .map(|&l| {
                // k=Begin -> l
                let p_init = self.node(l).init_prob();
                // emission prob on l
                let p_emit = self.p_match_emit(l, emission);
                p_init * ((param.p_MM * p_emit * t1.m[l]) + (param.p_MD * t0.d[l]))
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
    fn bib(&self, t0: &mut PHMMTable, t1: &PHMMTable, emission: u8, nodes: &[NodeIndex]) {
        let param = &self.param;

        // (1) to match and del of all nodes
        let p_to_match_del: Prob = nodes
            .iter()
            .map(|&l| {
                // k=Begin -> l
                let p_init = self.node(l).init_prob();
                // emission prob on l
                let p_emit = self.p_match_emit(l, emission);
                p_init * ((param.p_IM * p_emit * t1.m[l]) + (param.p_ID * t0.d[l]))
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
    fn be(&self, t0: &mut PHMMTable, _t1: &PHMMTable, _emission: u8, _nodes: &[NodeIndex]) {
        t0.e = Prob::from_prob(0.0);
    }
}

#[cfg(test)]
mod tests {}
