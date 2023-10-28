//!
//! Forward algorithm definitions
//!

use super::common::{PHMMEdge, PHMMModel, PHMMNode};
use super::hint::Mapping;
use super::table::{PHMMKind, PHMMTable, PHMMTables};
use crate::common::collection::Bases;
use crate::prob::{p, Prob};
use petgraph::graph::NodeIndex;

///
/// Forward Algorithm
///
impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    /// Run Forward algorithm to the emissions Dense
    ///
    /// `ft_i[k]` = P(emits `x[:i+1] = x[0],...,x[i]` and now in state `t_k`)
    ///
    /// * `t` is a type of state, either Match, Ins, Del
    /// * `k` is a node index
    ///
    pub fn forward<X: AsRef<Bases>>(&self, emissions: X) -> PHMMTables {
        let r0 = PHMMTables {
            init_table: self.f_init(true),
            tables: Vec::new(),
            kind: PHMMKind::Forward,
        };
        let all_nodes = self.to_all_nodes();
        emissions
            .as_ref()
            .into_iter()
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
    /// Run Forward algorithm to the emissions using hint/mapping information
    ///
    /// To calculate F.tables[i] = F[i+1], use S[i]=F[i+1]B[i+1] (0<=i<n)
    ///
    pub fn forward_with_mapping<X: AsRef<Bases>>(
        &self,
        emissions: X,
        mapping: &Mapping,
    ) -> PHMMTables {
        let r0 = PHMMTables {
            init_table: self.f_init(true),
            tables: Vec::new(),
            kind: PHMMKind::Forward,
        };
        emissions
            .as_ref()
            .into_iter()
            .enumerate()
            .fold(r0, |mut r, (i, &emission)| {
                let table_prev = if i == 0 {
                    &r.init_table
                } else {
                    r.last_table()
                };
                let table = self.f_step(i, emission, table_prev, mapping.nodes(i), false, false);
                r.tables.push(table);
                r
            })
    }
    ///
    /// Run Forward algorithm to the emissions using hint/mapping information and returns score only (drops intermediate result PHMMTable)
    ///
    pub fn forward_with_mapping_score_only<X: AsRef<Bases>>(
        &self,
        emissions: X,
        mapping: &Mapping,
    ) -> Prob {
        let mut table = self.f_init(true);
        for (i, &emission) in emissions.as_ref().into_iter().enumerate() {
            table = self.f_step(i, emission, &table, mapping.nodes(i), false, false);
        }
        table.e
    }
    ///
    /// Run Forward algorithm to the emissions, with sparse calculation
    ///
    pub fn forward_sparse<X: AsRef<Bases>>(&self, emissions: X, use_max_ratio: bool) -> PHMMTables {
        let r0 = PHMMTables {
            init_table: self.f_init(true),
            tables: Vec::new(),
            kind: PHMMKind::Forward,
        };
        let param = &self.param;
        let all_nodes = self.to_all_nodes();
        emissions
            .as_ref()
            .into_iter()
            .enumerate()
            .fold(r0, |mut r, (i, &emission)| {
                // previous table and its top_nodes
                let table_prev = if i == 0 {
                    &r.init_table
                } else {
                    r.last_table()
                };
                let top_nodes = if use_max_ratio {
                    table_prev.top_nodes_by_score_ratio(param.active_node_max_ratio)
                } else {
                    table_prev.top_nodes(param.n_active_nodes)
                };

                // the resulting table will be sparse/dense?
                let use_dense = if use_max_ratio {
                    if table_prev.is_dense() {
                        // if active_nodes is too small, switch to sparse.
                        // otherwise, use dense.
                        if i == 0 {
                            true
                        } else if i < param.n_warmup {
                            top_nodes.len() > param.warmup_threshold
                        } else {
                            false
                        }
                    } else {
                        // if previous table is already sparse, use sparse too.
                        false
                    }
                } else {
                    // not adaptive
                    // first n_warmup tables will be computed using dense.
                    i < param.n_warmup
                };

                // println!("i={} L={} use_dense={}", i, top_nodes.len(), use_dense);

                let table = if use_dense {
                    // dense_table
                    self.f_step(i, emission, table_prev, &all_nodes, true, false)
                } else {
                    // sparse_table
                    // determine next active nodes
                    let active_nodes = self.to_childs_and_us(&top_nodes);
                    self.f_step(i, emission, table_prev, &active_nodes, false, true)
                };
                r.tables.push(table);
                r
            })
    }
    ///
    ///
    ///
    pub fn forward_sparse_score_only<X: AsRef<Bases>>(
        &self,
        emissions: X,
        use_max_ratio: bool,
    ) -> Prob {
        let mut table = self.f_init(true);
        let param = &self.param;
        let all_nodes = self.to_all_nodes();
        for (i, &emission) in emissions.as_ref().into_iter().enumerate() {
            // previous top_nodes
            let top_nodes = if use_max_ratio {
                table.top_nodes_by_score_ratio(param.active_node_max_ratio)
            } else {
                table.top_nodes(param.n_active_nodes)
            };

            // the resulting table will be sparse/dense?
            let use_dense = if use_max_ratio {
                if table.is_dense() {
                    // if active_nodes is too small, switch to sparse.
                    // otherwise, use dense.
                    if i == 0 {
                        true
                    } else if i < param.n_warmup {
                        top_nodes.len() > param.warmup_threshold
                    } else {
                        false
                    }
                } else {
                    // if previous table is already sparse, use sparse too.
                    false
                }
            } else {
                // not adaptive
                // first n_warmup tables will be computed using dense.
                i < param.n_warmup
            };

            table = if use_dense {
                // dense_table
                self.f_step(i, emission, &table, &all_nodes, true, false)
            } else {
                // sparse_table
                let active_nodes = self.to_childs_and_us(&top_nodes);
                self.f_step(i, emission, &table, &active_nodes, false, true)
            };
        }
        table.e
    }
    ///
    /// Run Forward algorithm to the emissions, with sparse calculation
    ///
    pub fn forward_sparse_v0<X: AsRef<Bases>>(
        &self,
        emissions: X,
        use_max_ratio: bool,
    ) -> PHMMTables {
        let r0 = PHMMTables {
            init_table: self.f_init(true),
            tables: Vec::new(),
            kind: PHMMKind::Forward,
        };
        let param = &self.param;
        let all_nodes = self.to_all_nodes();
        emissions
            .as_ref()
            .into_iter()
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
                    let active_nodes = if use_max_ratio {
                        self.to_childs_and_us(
                            &table_prev.top_nodes_by_score_ratio(param.active_node_max_ratio),
                        )
                    } else {
                        self.to_childs_and_us(&table_prev.top_nodes(param.n_active_nodes))
                    };
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
    /// * is_dense
    ///     if true, the resulting PHMMTable will be dense. (Store probabilities for all nodes)
    /// * is_adaptive
    ///     use adaptive node calculation in Del states
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
mod tests {
    use super::*;
    use crate::common::{ni, sequence_to_string};
    use crate::e2e;
    use crate::hmmv2::mocks::*;
    use crate::hmmv2::params::PHMMParams;
    use crate::kmer::VecKmer;
    use crate::prob::lp;
    use approx::{assert_abs_diff_eq, assert_relative_eq};
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
        let r2 = phmm.forward_sparse(&read, false);
        for i in 0..r1.n_emissions() {
            let t1 = r1.table(i);
            let t2 = r2.table(i);
            println!("i={}", i);
            println!("{}", t1);
            println!("{}", t2);
            let d = t1.diff(&t2);
            println!("{}", d);
            assert!(d < 0.000000001);
        }
    }
    #[test]
    fn hmm_forward_with_hint_mock_linear_high_error() {
        let phmm = mock_linear_phmm(PHMMParams::high_error());
        // read1
        let read1 = b"CGATC";
        let o = phmm.run(read1);
        let hint = o.to_mapping(3);
        println!("{:?}", hint);
        println!("{:?}", hint.len());
        println!("{}", hint);
        assert_eq!(
            hint.nodes,
            vec![
                vec![ni(3), ni(2), ni(4)],
                vec![ni(4), ni(3), ni(5)],
                vec![ni(5), ni(6), ni(4)],
                vec![ni(6), ni(7), ni(5)],
                vec![ni(7), ni(8), ni(6)],
            ]
        );

        let r1 = phmm.forward(read1);
        let r2 = phmm.forward_with_mapping(read1, &hint);
        println!("{}", r1.last_table());
        println!("{}", r2.last_table());
        let p1 = r1.full_prob();
        let p2 = r2.full_prob();
        println!("p(dense)={}", p1);
        println!("p(hint)={}", p2);
        assert!(p1.log_diff(p2) < 0.1);
    }
}
