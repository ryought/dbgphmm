use super::common::PHMMModel;
use super::table::{PHMMTable, PHMMTables};
use crate::prob::{p, Prob};
use petgraph::graph::NodeIndex;

///
/// Forward algorithm
///
impl PHMMModel {
    ///
    /// Run Forward algorithm to the emissions Dense
    ///
    /// `ft_i[k]` = P(emits `x[:i+1] = x[0],...,x[i]` and now in state `t_k`)
    ///
    /// * `t` is a type of state, either Match, Ins, Del
    /// * `k` is a node index
    ///
    pub fn forward(&self, emissions: &[u8]) -> PHMMTables {
        let r0 = PHMMTables {
            init_table: self.f_init(true),
            tables: Vec::new(),
            is_forward: true,
        };
        let all_nodes = self.to_all_nodes();
        emissions
            .iter()
            .enumerate()
            .fold(r0, |mut r, (i, &emission)| {
                let table_prev = if i == 0 {
                    &r.init_table
                } else {
                    r.last_table()
                };
                let table = self.f_step(i, emission, table_prev, &all_nodes, true, false);
                r.tables.push(table);
                r
            })
    }
    ///
    /// Run Forward algorithm to the emissions, with sparse calculation
    ///
    pub fn forward_sparse(&self, emissions: &[u8]) -> PHMMTables {
        let r0 = PHMMTables {
            init_table: self.f_init(true),
            tables: Vec::new(),
            is_forward: true,
        };
        let param = &self.param;
        let all_nodes = self.to_all_nodes();
        emissions
            .iter()
            .enumerate()
            .fold(r0, |mut r, (i, &emission)| {
                if i < param.n_warmup {
                    // dense_table
                    let table_prev = if i == 0 {
                        &r.init_table
                    } else {
                        r.last_table()
                    };
                    let table = self.f_step(i, emission, table_prev, &all_nodes, true, false);
                    r.tables.push(table);
                } else {
                    // sparse_table
                    let table_prev = r.last_table();
                    let active_nodes = self.to_childs(&table_prev.top_nodes(param.n_active_nodes));
                    let table = self.f_step(i, emission, table_prev, &active_nodes, false, true);
                    r.tables.push(table);
                };
                r
            })
    }
    ///
    /// Create init_table in PHMMResult for Forward algorithm
    ///
    fn f_init(&self, is_dense: bool) -> PHMMTable {
        PHMMTable::new(
            is_dense,
            self.n_nodes(),
            p(0.0),
            p(0.0),
            p(0.0),
            p(1.0), // only MatchBegin has probability (p=1)
            p(0.0),
            p(0.0),
        )
    }
    ///
    /// Calculate the table from the previous table
    /// for Forward algorithm
    ///
    fn f_step(
        &self,
        _i: usize,
        emission: u8,
        prev_table: &PHMMTable,
        active_nodes: &[NodeIndex],
        is_dense: bool,
        is_adaptive: bool,
    ) -> PHMMTable {
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

        // normal state first
        self.fm(&mut table, prev_table, emission, active_nodes);
        self.fi(&mut table, prev_table, emission, active_nodes);
        self.fmb(&mut table, prev_table, emission, active_nodes);
        self.fib(&mut table, prev_table, emission, active_nodes);
        // silent state next
        self.fd(&mut table, prev_table, emission, active_nodes, is_adaptive);
        self.fe(&mut table, prev_table, emission, active_nodes);

        table
    }
}

// functions to calculate each step
impl PHMMModel {
    /// Fill the forward probs of `Match` states
    ///
    /// For `i=0,...,n-1`
    ///
    /// ```text
    /// fm_i[k]
    /// = P(emits x[:i+1] and in state m_k)
    /// =   (from_m_parents) \sum_{l: parents} t_lk p_mm e(x[i]) fm_i-1[l]
    ///   + (from_i_parents) \sum_{l: parents} t_lk p_im e(x[i]) fi_i-1[l]
    ///   + (from_d_parents) \sum_{l: parents} t_lk p_dm e(x[i]) fd_i-1[l]
    ///   + (from_m_begin)                     t_bk p_mm e(x[i]) fm_i-1[b]
    ///   + (from_i_begin)                     t_bk p_im e(x[i]) fi_i-1[b]
    ///
    /// = e(x[i]) * (
    ///       \sum_{l: parents} t_lk (p_mm fm_i-1[l] + p_im fi_i-1[l] + p_dm fd_i-1[l])
    ///       +
    ///       t_bk (p_mm fm_i-1[b] + p_im fi_i-1[b])
    ///   )
    /// ```
    ///
    /// (Here `x[:i+1] = x[0],...,x[i]`)
    ///
    /// calculate `t0.m` from `t1.m, t1.i, t1.d, t1.mb, t1.ib`
    fn fm(&self, t0: &mut PHMMTable, t1: &PHMMTable, emission: u8, nodes: &[NodeIndex]) {
        let param = &self.param;
        for &k in nodes {
            // emission prob
            let p_emit = self.p_match_emit(k, emission);
            // init_prob
            let p_init = self.node(k).init_prob();

            // (1) from normal node
            let from_normal: Prob = self
                .parents(k)
                .map(|(_, l, ew)| {
                    let p_trans = ew.trans_prob();
                    p_trans * (param.p_MM * t1.m[l] + param.p_IM * t1.i[l] + param.p_DM * t1.d[l])
                })
                .sum();

            // (2) from begin node
            let from_begin = p_init * (param.p_MM * t1.mb + param.p_IM * t1.ib);

            t0.m[k] = p_emit * (from_normal + from_begin);
        }
    }
    /// Fill the forward probs of `Ins` states
    ///
    /// For `i=0,...,n-1`
    ///
    /// ```text
    /// fm_i[k]
    /// = P(emits x[:i+1] and then in state i_k)
    /// =   (from_self_m) p_mi e(x[i]) fm_i-1[k]
    ///   + (from_self_i) p_ii e(x[i]) fi_i-1[k]
    ///   + (from_self_d) p_di e(x[i]) fd_i-1[k]
    /// = e(x[i]) * (
    ///         p_mi fm_i-1[k] + p_ii fi_i-1[k] + p_di fd_i-1[k]
    ///   )
    /// ```
    ///
    /// (Here `x[:i+1] = x[0],...,x[i]`)
    ///
    /// calculate `t0.i` from `t1.m, t1.i, t1.d, t1.mb, t1.ib`
    fn fi(&self, t0: &mut PHMMTable, t1: &PHMMTable, _emission: u8, nodes: &[NodeIndex]) {
        let param = &self.param;
        for &k in nodes {
            // emission prob
            let p_emit = self.p_ins_emit();

            // from my own node
            let from_me = param.p_MI * t1.m[k] + param.p_II * t1.i[k] + param.p_DI * t1.d[k];
            t0.i[k] = p_emit * from_me;
        }
    }
    /// fill the forward probs of `Del` states
    ///
    /// For `i=0,...,n-1`
    ///
    /// ```text
    /// fd_i[k]
    /// = P(emits x[:i+1] and then in state d_k)
    /// =   (from_m_parents) \sum_{l: parents} t_lk p_md fm_i[l]
    ///   + (from_i_parents) \sum_{l: parents} t_lk p_id fi_i[l]
    ///   + (from_d_parents) \sum_{l: parents} t_lk p_dd fd_i[l]
    ///   + (from_m_begin)                     t_bk p_md fm_i[b]
    ///   + (from_i_begin)                     t_bk p_id fi_i[b]
    /// ```
    ///
    /// (Here `x[:i+1] = x[0],...,x[i]`)
    ///
    /// This calculation has a recursive definition of `fd_i`, so purging them
    /// by allowing only `param.n_max_gaps` times continuous deletions.
    ///
    /// ```text
    /// fd_i(0)[k]
    /// =   (from_m_parents) \sum_{l: parents} t_lk p_md fm_i[l]
    ///   + (from_i_parents) \sum_{l: parents} t_lk p_id fi_i[l]
    ///   + (from_m_begin)                     t_bk p_md fm_i[b]
    ///   + (from_i_begin)                     t_bk p_id fi_i[b]
    ///
    /// t = 1,2,3,...
    /// fd_i(t)[k]
    /// =   (from_d_parents) \sum_{l: parents} t_lk p_dd fd_i(t-1)[l]
    ///
    /// fd_i[k] = \sum_t fd_i(t)[k]
    /// ```
    ///
    /// calculate t0.d from t0.m and t0.i
    fn fd(
        &self,
        t0: &mut PHMMTable,
        _t1: &PHMMTable,
        _emission: u8,
        nodes: &[NodeIndex],
        is_adaptive: bool,
    ) {
        let param = &self.param;
        let mut active_nodes = None;
        let is_dense = t0.is_dense();

        // run t=0
        if is_adaptive {
            active_nodes = Some(self.to_childs(&nodes));
        }
        let mut fdt0 = self.fd0(
            t0,
            if is_adaptive {
                active_nodes.as_ref().unwrap()
            } else {
                nodes
            },
            is_dense,
        );
        t0.d += &fdt0.d;

        for _t in 0..param.n_max_gaps {
            // run t+1
            if is_adaptive {
                active_nodes = Some(self.to_childs(active_nodes.as_ref().unwrap()));
            }
            fdt0 = self.fdt(
                &fdt0,
                if is_adaptive {
                    active_nodes.as_ref().unwrap()
                } else {
                    nodes
                },
                is_dense,
            );
            t0.d += &fdt0.d;
        }
    }
    ///
    /// Calculate `fd_i(t=0)[k]`
    ///
    /// ```text
    /// fd_i(0)[k]
    /// =   (from_m_parents) \sum_{l: parents} t_lk p_md fm_i[l]
    ///   + (from_i_parents) \sum_{l: parents} t_lk p_id fi_i[l]
    ///   + (from_m_begin)                     t_bk p_md fm_i[b]
    ///   + (from_i_begin)                     t_bk p_id fi_i[b]
    /// ```
    ///
    /// active_nodes will be determined by the childs of t0
    ///
    fn fd0(&self, t0: &PHMMTable, nodes: &[NodeIndex], is_dense: bool) -> PHMMTable {
        let param = &self.param;
        let mut fd0 = PHMMTable::zero(is_dense, self.n_nodes());
        for &k in nodes {
            // (1) from normal node
            let from_normal: Prob = self
                .parents(k)
                .map(|(_, l, ew)| {
                    // l -> k
                    let p_trans = ew.trans_prob();
                    p_trans * (param.p_MD * t0.m[l] + param.p_ID * t0.i[l])
                })
                .sum();

            // (2) from begin node
            let p_init = self.node(k).init_prob();
            let from_begin = p_init * (param.p_MD * t0.mb + param.p_ID * t0.ib);

            fd0.d[k] = from_normal + from_begin;
        }
        fd0
    }
    ///
    /// Calculate `fd_i(t)[k]` for `t > 0`
    ///
    /// ```text
    /// t = 1,2,3,...
    /// fd_i(t)[k]
    /// =   (from_d_parents) \sum_{l: parents} t_lk p_dd fd_i(t-1)[l]
    /// ```
    fn fdt(&self, fdt1: &PHMMTable, nodes: &[NodeIndex], is_dense: bool) -> PHMMTable {
        let param = &self.param;
        let mut fdt0 = PHMMTable::zero(is_dense, self.n_nodes());
        for &k in nodes {
            fdt0.d[k] = self
                .parents(k)
                .map(|(_, l, ew)| {
                    // l -> k
                    let p_trans = ew.trans_prob();
                    p_trans * (param.p_DD * fdt1.d[l])
                })
                .sum();
        }
        fdt0
    }
    /// fill `MatchBegin` states
    ///
    /// ```text
    /// fm_i[b] = 1 (if i==-1)
    ///           0 (otherwise)
    /// ```
    fn fmb(&self, t0: &mut PHMMTable, _t1: &PHMMTable, _emission: u8, _nodes: &[NodeIndex]) {
        t0.mb = Prob::from_prob(0.0);
    }
    /// fill `InsBegin` states
    ///
    /// ```text
    /// fi_i[b]
    /// =   (from_self_m) p_mi e(x[i]) fm_i-1[b]
    ///   + (from_self_i) p_ii e(x[i]) fi_i-1[b]
    /// ```
    fn fib(&self, t0: &mut PHMMTable, t1: &PHMMTable, _emission: u8, _nodes: &[NodeIndex]) {
        let param = &self.param;
        let p_emit = self.p_ins_emit();
        t0.ib = p_emit * (param.p_MI * t1.mb + param.p_II * t1.ib);
    }
    /// fill `End` state
    ///
    /// ```text
    /// fe_i
    /// =   (from_all_m) \sum_{k} p_e fm_i[k]
    ///   + (from_all_i) \sum_{k} p_e fi_i[k]
    ///   + (from_all_d) \sum_{k} p_e fd_i[k]
    /// ```
    fn fe(&self, t0: &mut PHMMTable, _t1: &PHMMTable, _emission: u8, nodes: &[NodeIndex]) {
        let param = &self.param;
        let p_normal: Prob = nodes.iter().map(|&k| t0.m[k] + t0.i[k] + t0.d[k]).sum();
        t0.e = param.p_end * p_normal
    }
}

//
// Tests
//

#[cfg(test)]
mod tests {}
