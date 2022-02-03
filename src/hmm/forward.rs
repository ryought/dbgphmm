use super::common::{PHMMEdge, PHMMModel, PHMMNode};
use super::table::{PHMMResult, PHMMTable};
use crate::prob::Prob;
use crate::veclike::{DenseVec, VecLike};

impl<N: PHMMNode, E: PHMMEdge> PHMMModel<N, E> {
    ///
    /// Run Forward algorithm to the emissions
    ///
    pub fn forward<V: VecLike<Prob>>(&self, emissions: &[u8]) -> PHMMResult<V> {
        let tables_init = vec![self.f_init()];
        let tables =
            emissions
                .iter()
                .enumerate()
                .fold(tables_init, |mut tables, (i, &emission)| {
                    let table = self.f_step(i, emission, tables.last().unwrap());
                    tables.push(table);
                    tables
                });
        PHMMResult(tables)
    }
    /// fill Del states
    /// calculate t0.d from t0.m and t0.i
    fn fd<V: VecLike<Prob>>(&self, t0: &mut PHMMTable<V>) {}
    fn f_init<V: VecLike<Prob>>(&self) -> PHMMTable<V> {
        let n_kmers = 10;
        PHMMTable {
            // only MatchBegin has probability (p=1)
            mb: Prob::from_prob(1.0),
            // the other states has p=0
            m: V::new(n_kmers, Prob::from_prob(0.0)),
            i: V::new(n_kmers, Prob::from_prob(0.0)),
            d: V::new(n_kmers, Prob::from_prob(0.0)),
            ib: Prob::from_prob(0.0),
            e: Prob::from_prob(0.0),
        }
    }
    fn f_step<V: VecLike<Prob>>(
        &self,
        i: usize,
        emission: u8,
        prev: &PHMMTable<V>,
    ) -> PHMMTable<V> {
        unimplemented!();
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
