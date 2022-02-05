//!
//! Forward algorithm definitions
//!

use super::common::{PHMMEdge, PHMMModel, PHMMNode};
use super::table::{PHMMResult, PHMMTable};
use crate::prob::{p, Prob};
use crate::vector::{NodeVec, Storage};

// wrappers and exposed functions
impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    /// Run Forward algorithm to the emissions
    ///
    pub fn forward<S>(&self, emissions: &[u8]) -> PHMMResult<S>
    where
        S: Storage<Item = Prob>,
    {
        let r0 = PHMMResult {
            init_table: self.f_init(),
            tables: Vec::new(),
        };
        emissions
            .iter()
            .enumerate()
            .fold(r0, |mut r, (i, &emission)| {
                let table = if i == 0 {
                    self.f_step(i, emission, &r.init_table)
                } else {
                    self.f_step(i, emission, r.tables.last().unwrap())
                };
                r.tables.push(table);
                r
            })
    }
    ///
    /// Create init_table in PHMMResult for Forward algorithm
    ///
    fn f_init<S>(&self) -> PHMMTable<S>
    where
        S: Storage<Item = Prob>,
    {
        PHMMTable::new(
            self.n_nodes(),
            p(0.0),
            p(0.0),
            p(0.0),
            // only MatchBegin has probability (p=1)
            p(1.0),
            p(0.0),
            p(0.0),
        )
    }
    ///
    /// Calculate the table from the previous table
    /// for Forward algorithm
    ///
    fn f_step<S>(&self, _i: usize, emission: u8, prev_table: &PHMMTable<S>) -> PHMMTable<S>
    where
        S: Storage<Item = Prob>,
    {
        let mut table = PHMMTable::new(
            self.n_nodes(),
            p(0.0),
            p(0.0),
            p(0.0),
            p(0.0),
            p(0.0),
            p(0.0),
        );
        // normal state first
        self.fm(&mut table, prev_table, emission);
        self.fi(&mut table, prev_table, emission);
        self.fmb(&mut table, prev_table, emission);
        self.fib(&mut table, prev_table, emission);
        // silent state next
        self.fd(&mut table, prev_table, emission);
        self.fe(&mut table, prev_table, emission);
        table
    }
}

// functions to calculate each step
impl<'a, N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
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
    fn fm<S>(&self, t0: &mut PHMMTable<S>, t1: &PHMMTable<S>, emission: u8)
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        for (k, kw) in self.nodes() {
            // emission prob
            let p_emit = if kw.emission() == emission {
                param.p_match
            } else {
                param.p_mismatch
            };

            // (1) from normal node
            let from_normal: Prob = self
                .parents(k)
                .map(|(_, l, ew)| {
                    let p_trans = ew.trans_prob();
                    p_trans * (param.p_MM * t1.m[l] + param.p_IM * t1.i[l] + param.p_DM * t1.d[l])
                })
                .sum();

            // (2) from begin node
            let from_begin = kw.init_prob() * (param.p_MM * t1.mb + param.p_IM * t1.ib);

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
    fn fi<S>(&self, t0: &mut PHMMTable<S>, t1: &PHMMTable<S>, _emission: u8)
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        for (k, _) in self.nodes() {
            // emission prob
            let p_emit = param.p_random;

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
    fn fd<S>(&self, t0: &mut PHMMTable<S>, _t1: &PHMMTable<S>, _emission: u8)
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        let mut fdt0 = self.fd0(t0);
        t0.d += &fdt0;
        for _t in 0..param.n_max_gaps {
            fdt0 = self.fdt(&fdt0);
            t0.d += &fdt0;
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
    fn fd0<S>(&self, t0: &PHMMTable<S>) -> NodeVec<S>
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        let mut fd0 = NodeVec::new(self.n_nodes(), Prob::from_prob(0.0));
        for (k, kw) in self.nodes() {
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
            let from_begin = kw.init_prob() * (param.p_MD * t0.mb + param.p_ID * t0.ib);

            fd0[k] = from_normal + from_begin;
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
    fn fdt<S>(&self, fdt1: &NodeVec<S>) -> NodeVec<S>
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        let mut fdt0 = NodeVec::new(self.n_nodes(), Prob::from_prob(0.0));
        for (k, _) in self.nodes() {
            fdt0[k] = self
                .parents(k)
                .map(|(_, l, ew)| {
                    // l -> k
                    let p_trans = ew.trans_prob();
                    p_trans * (param.p_DD * fdt1[l])
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
    fn fmb<S>(&self, t0: &mut PHMMTable<S>, _t1: &PHMMTable<S>, _emission: u8)
    where
        S: Storage<Item = Prob>,
    {
        t0.mb = Prob::from_prob(0.0);
    }
    /// fill `InsBegin` states
    ///
    /// ```text
    /// fi_i[b]
    /// =   (from_self_m) p_mi e(x[i]) fm_i-1[b]
    ///   + (from_self_i) p_ii e(x[i]) fi_i-1[b]
    /// ```
    fn fib<S>(&self, t0: &mut PHMMTable<S>, t1: &PHMMTable<S>, _emission: u8)
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        let p_emit = param.p_random;
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
    fn fe<S>(&self, t0: &mut PHMMTable<S>, _t1: &PHMMTable<S>, _emission: u8)
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        let p_normal: Prob = self.nodes().map(|(k, _)| t0.m[k] + t0.i[k] + t0.d[k]).sum();
        t0.e = param.p_end * p_normal
    }
}

#[cfg(test)]
mod tests {
    use super::super::seqgraph::create_linear_seq_graph;

    #[test]
    fn create_linear_seq_graph_test() {
        let g = create_linear_seq_graph(b"ATCGGCTAGC");
        let phmm = g.to_phmm();
        println!("{}", phmm);
    }
}
