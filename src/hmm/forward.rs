use super::common::{PHMMEdge, PHMMModel, PHMMNode};
use super::table::{PHMMResult, PHMMTable};
use super::veclikewrap::NodeVec;
use crate::prob::Prob;
use crate::veclike::VecLike;

// wrappers and exposed functions
impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    /// Run Forward algorithm to the emissions
    ///
    pub fn forward<V: VecLike<Prob>>(&self, emissions: &[u8]) -> PHMMResult<V> {
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
    fn f_init<V: VecLike<Prob>>(&self) -> PHMMTable<V> {
        PHMMTable::new(
            self.n_nodes(),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
            // only MatchBegin has probability (p=1)
            Prob::from_prob(1.0),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
        )
    }
    ///
    /// Calculate the table from the previous table
    /// for Forward algorithm
    ///
    fn f_step<V: VecLike<Prob>>(
        &self,
        i: usize,
        emission: u8,
        prev_table: &PHMMTable<V>,
    ) -> PHMMTable<V> {
        let mut table = PHMMTable::new(
            self.n_nodes(),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
            Prob::from_prob(0.0),
        );
        self.fm(&mut table, prev_table, emission);
        self.fi(&mut table, prev_table, emission);
        self.fb(&mut table, prev_table, emission);
        self.fd(&mut table, prev_table, emission);
        self.fe(&mut table, prev_table, emission);
        table
    }
}

// functions to calculate each step
impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    /// fill `Match` states
    /// calculate `t0.m` from `t1.m, t1.i, t1.d`
    fn fm<V: VecLike<Prob>>(&self, t0: &mut PHMMTable<V>, t1: &PHMMTable<V>, emission: u8) {
        let param = &self.param;
        for (k, kw) in self.nodes() {
            // emission prob
            let p_emit = if kw.emission() == emission {
                param.p_match
            } else {
                param.p_mismatch
            };

            // from normal node
            let p_normal: Prob = self
                .parents(k)
                .map(|(_, l, ew)| {
                    let p_trans = ew.trans_prob();
                    p_trans * (param.p_MM * t1.m[l] + param.p_IM * t1.i[l] + param.p_DM * t1.d[l])
                })
                .sum();
            // from begin node
            let p_begin = kw.init_prob() * (param.p_MM * t1.mb + param.p_IM * t1.ib);
            t0.m[k] = p_emit * (p_normal + p_begin);
        }
    }
    /// fill `Ins` states
    /// calculate `t0.i` from `t1.m, t1.i, t1.d`
    fn fi<V: VecLike<Prob>>(&self, t0: &mut PHMMTable<V>, t1: &PHMMTable<V>, emission: u8) {
        let param = &self.param;
        for (k, kw) in self.nodes() {
            // emission prob
            let p_emit = param.p_random;

            // from my own node
            let p_normal = param.p_MI * t1.m[k] + param.p_II * t1.i[k] + param.p_DI * t1.d[k];
            t0.i[k] = p_emit * p_normal;
        }
    }
    /// fill `Del` states
    /// calculate t0.d from t0.m and t0.i
    fn fd<V: VecLike<Prob>>(&self, t0: &mut PHMMTable<V>, _t1: &PHMMTable<V>, _emission: u8) {
        let param = &self.param;
        // TODO
        let mut fdt0 = self.fd0(t0);
        for t in 0..param.n_max_gaps {
            let fdt1 = self.fdt(&fdt0);
        }
    }
    fn fd0<V: VecLike<Prob>>(&self, t0: &PHMMTable<V>) -> NodeVec<V> {
        let param = &self.param;
        let mut fd0 = NodeVec::new(self.n_nodes(), Prob::from_prob(0.0));
        for (k, kw) in self.nodes() {
            let p_normal: Prob = self
                .parents(k)
                .map(|(_, l, ew)| {
                    let p_trans = ew.trans_prob();
                    p_trans * (param.p_MD * t0.m[l] + param.p_ID * t0.i[l])
                })
                .sum();
            let p_begin = kw.init_prob() * (param.p_MD * t0.m[k] + param.p_ID * t0.i[k]);
            fd0[k] = p_normal + p_begin;
        }
        fd0
    }
    fn fdt<V: VecLike<Prob>>(&self, fdt1: &NodeVec<V>) -> NodeVec<V> {
        let param = &self.param;
        let mut fdt0 = NodeVec::new(self.n_nodes(), Prob::from_prob(0.0));
        for (k, kw) in self.nodes() {
            fdt0[k] = self
                .parents(k)
                .map(|(_, l, ew)| {
                    let p_trans = ew.trans_prob();
                    p_trans * (param.p_DD * fdt1[l])
                })
                .sum();
        }
        fdt0
    }
    /// fill `MatchBegin` and `InsBegin` states
    fn fb<V: VecLike<Prob>>(&self, t0: &mut PHMMTable<V>, t1: &PHMMTable<V>, _emission: u8) {
        let param = &self.param;
        let p_emit = param.p_random;
        t0.mb = Prob::from_prob(0.0);
        t0.ib = p_emit * (param.p_MI * t1.mb + param.p_II * t1.ib);
    }
    /// fill `End` state
    fn fe<V: VecLike<Prob>>(&self, t0: &mut PHMMTable<V>, t1: &PHMMTable<V>, _emission: u8) {
        let param = &self.param;
        let p_normal: Prob = self.nodes().map(|(k, _)| t1.m[k] + t1.i[k] + t1.d[k]).sum();
        t0.e = param.p_end * p_normal
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
