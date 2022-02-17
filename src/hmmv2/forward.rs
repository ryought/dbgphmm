//!
//! Forward algorithm definitions
//!

use super::common::{PHMMEdge, PHMMModel, PHMMNode};
use super::result::{PHMMResult, PHMMResultLike, PHMMResultSparse};
use super::table::PHMMTable;
use super::table_ref::PHMMTableRef;
use crate::prob::{p, Prob};
use crate::vector::{NodeVec, Storage};

// wrappers and exposed functions
impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    /// Run Forward algorithm to the emissions
    ///
    /// `ft_i[k]` = P(emits `x[:i+1] = x[0],...,x[i]` and now in state `t_k`)
    ///
    /// * `t` is a type of state, either Match, Ins, Del
    /// * `k` is a node index
    ///
    pub fn forward(&self, emissions: &[u8]) -> PHMMResult {
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
    /// Run Forward algorithm to the emissions, with sparse calculation
    /// Not tested
    ///
    pub fn forward_sparse(
        &self,
        emissions: &[u8],
        n_warmup: usize,
        n_active_nodes: usize,
    ) -> PHMMResultSparse {
        assert!(n_warmup > 1);
        assert!(n_active_nodes > 1);
        let r0 = PHMMResultSparse {
            init_table: self.f_init(),
            tables_warmup: Vec::new(),
            tables_sparse: Vec::new(),
            is_forward: true,
        };
        emissions
            .iter()
            .enumerate()
            .fold(r0, |mut r, (i, &emission)| {
                if i < n_warmup {
                    // dense_table -> dense_table
                    let table = if i == 0 {
                        self.f_step(i, emission, &r.init_table)
                    } else {
                        self.f_step(i, emission, r.tables_warmup.last().unwrap())
                    };
                    r.tables_warmup.push(table);
                } else if i == n_warmup {
                    // dense_table -> sparse_table
                    let table_prev = r
                        .tables_warmup
                        .last()
                        .unwrap()
                        .to_sparse_active_nodes(n_active_nodes);
                    let mut table = self.f_step(i, emission, &table_prev);
                    table.refresh_active_nodes(n_active_nodes);
                    r.tables_sparse.push(table);
                } else {
                    // sparse_table -> sparse_table
                    let mut table = self.f_step(i, emission, r.tables_sparse.last().unwrap());
                    table.refresh_active_nodes(n_active_nodes);
                    r.tables_sparse.push(table);
                };
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
    /// If `prev_table.active_nodes` is set (i.e. `ActiveNodes::Only`),
    ///
    fn f_step<S>(&self, _i: usize, emission: u8, prev_table: &PHMMTable<S>) -> PHMMTable<S>
    where
        S: Storage<Item = Prob>,
    {
        // candidates of active nodes of next step
        let active_nodes = prev_table.active_nodes.to_childs(self);

        let mut table = PHMMTable::new_with_active_nodes(
            self.n_nodes(),
            p(0.0),
            p(0.0),
            p(0.0),
            p(0.0),
            p(0.0),
            p(0.0),
            active_nodes,
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
    /// If `t0.active_nodes` is set, only node k in the ActiveNodes will
    /// be calculated.
    ///
    /// calculate `t0.m` from `t1.m, t1.i, t1.d, t1.mb, t1.ib`
    fn fm<S>(&self, t0: &mut PHMMTable<S>, t1: &PHMMTable<S>, emission: u8)
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        for (k, kw) in self.active_nodes(&t0.active_nodes) {
            // emission prob
            let p_emit = self.p_match_emit(k, emission);

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
        for (k, _) in self.active_nodes(&t0.active_nodes) {
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
    fn fd<S>(&self, t0: &mut PHMMTable<S>, _t1: &PHMMTable<S>, _emission: u8)
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        let mut fdt0 = self.fd0(t0);
        t0.d += &fdt0.d;
        for _t in 0..param.n_max_gaps {
            fdt0 = self.fdt(&fdt0);
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
    fn fd0<S>(&self, t0: &PHMMTable<S>) -> PHMMTable<S>
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        let active_nodes = t0.active_nodes.to_childs(self);
        let mut fd0 = PHMMTable::zero_with_active_nodes(self.n_nodes(), active_nodes);
        for (k, kw) in self.active_nodes(&fd0.active_nodes) {
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
    fn fdt<S>(&self, fdt1: &PHMMTable<S>) -> PHMMTable<S>
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        let active_nodes = fdt1.active_nodes.to_childs(self);
        let mut fdt0 = PHMMTable::zero_with_active_nodes(self.n_nodes(), active_nodes);
        for (k, _) in self.active_nodes(&fdt0.active_nodes) {
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
    fn fe<S>(&self, t0: &mut PHMMTable<S>, _t1: &PHMMTable<S>, _emission: u8)
    where
        S: Storage<Item = Prob>,
    {
        let param = &self.param;
        let p_normal: Prob = self
            .active_nodes(&t0.active_nodes)
            .map(|(k, _)| t0.m[k] + t0.i[k] + t0.d[k])
            .sum();
        t0.e = param.p_end * p_normal
    }
}

//
// Tests
//

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{ni, sequence_to_string};
    use crate::hmm::params::PHMMParams;
    use crate::hmmv2::mocks::*;
    use crate::prob::lp;
    use crate::vector::DenseStorage;
    #[test]
    fn hmm_forward_mock_linear_zero_error() {
        let phmm = mock_linear_phmm(PHMMParams::zero_error());
        let r = phmm.forward(b"CGATC");
        assert_eq!(r.tables.len(), 5);
        assert_abs_diff_eq!(r.tables[2].m[ni(5)], lp(-2.3026250931), epsilon = 0.00001);
        assert_abs_diff_eq!(r.tables[3].m[ni(6)], lp(-2.3026250931), epsilon = 0.00001);
        assert_abs_diff_eq!(r.tables[4].m[ni(7)], lp(-2.3026350932), epsilon = 0.00001);
        // total probability
        assert_abs_diff_eq!(r.tables[4].e, lp(-13.8155605), epsilon = 0.00001);
        // no insertion and deletions, so i/d should be 0.
        for table in r.tables {
            for i in 0..table.n_nodes() {
                assert!(table.i[ni(i)].is_zero());
                assert!(table.d[ni(i)].is_zero());
            }
        }
        // with allowing no errors, CGATT cannot be emitted.
        // so it should have p=0
        let r2 = phmm.forward(b"CGATT");
        assert_eq!(r2.tables.len(), 5);
        assert!(r2.tables[4].e.is_zero());
    }
    #[test]
    fn hmm_forward_mock_linear_high_error() {
        let phmm = mock_linear_phmm(PHMMParams::high_error());
        // read 1
        let r = phmm.forward(b"CGATC");
        for table in r.tables.iter() {
            println!("{}", table);
        }
        println!("{}", r.init_table);
        assert_eq!(r.tables.len(), 5);
        assert_abs_diff_eq!(r.tables[4].e, lp(-15.212633254), epsilon = 0.00001);
        assert_abs_diff_eq!(r.tables[4].m[ni(7)], lp(-3.8652938682), epsilon = 0.00001);
        // read 2
        let r2 = phmm.forward(b"CGATT");
        assert_abs_diff_eq!(r2.tables[4].e, lp(-16.7862972), epsilon = 0.00001);
        // r[:4] and r2[:4] is the same emissions
        assert_abs_diff_eq!(r2.tables[3].e, r.tables[3].e, epsilon = 0.00001);
        assert_eq!(r2.tables.len(), 5);
        for table in r2.tables.iter() {
            println!("{}", table);
        }
    }
    #[test]
    fn hmm_forward_mock_sparse() {
        let phmm = mock_linear_random_phmm(100, 0, PHMMParams::default());
        let read = phmm.sample_read(32, 0);
        println!("{}", sequence_to_string(&read));

        let r1 = phmm.forward(&read);
        let r2 = phmm.forward_sparse(&read, 16, 40);
        for i in 0..r1.n_emissions() {
            let t1 = r1.table(i);
            let t2 = r2.table(i);
            let d = t1.diff(&t2);
            println!("{}", d);
            assert!(d < 0.000000001);
        }
    }
}
