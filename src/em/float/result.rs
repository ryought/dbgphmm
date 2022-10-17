//!
//! result structures
//!
use crate::dbg::float::{q_score_diff_exact, CopyDensity, FloatDbg, FloatDbgEdge, FloatDbgNode};
use crate::hmmv2::q::{q_score_exact, QScore};
use crate::prelude::*;
use crate::vector::{DenseStorage, NodeVec};

//
// for whole
//

pub type StepResult<K> = (
    (FloatDbg<K>, Prob), // init
    (FloatDbg<K>, Prob), // optimized
    (FloatDbg<K>, Prob), // shrinked
);

//
// For EM
//

pub struct EMResult<K: KmerLike> {
    pub e: Vec<Prob>,
    pub m: Vec<Vec<MStepResult<K>>>,
}

impl<K: KmerLike> EMResult<K> {
    pub fn new() -> Self {
        EMResult {
            e: Vec::new(),
            m: Vec::new(),
        }
    }
    pub fn to_final_dbg(&self) -> Option<FloatDbg<K>> {
        for (em_id, m_step_result) in self.m.iter().rev().enumerate() {
            for (m_id, m_step_once_result) in m_step_result.iter().rev().enumerate() {
                match m_step_once_result {
                    MStepResult::Init(dbg) => {
                        return Some(dbg.clone());
                    }
                    MStepResult::Update(dbg, _) => {
                        return Some(dbg.clone());
                    }
                    _ => {}
                };
            }
        }
        return None;
    }
    ///
    /// create Vec<NodeAttrVec> (that represents the time series of copy densities of nodes of EM
    /// steps) by using true dbg (Dbg<N, E>, not floated) and result (EMResult<K>).
    ///
    pub fn to_node_historys(&self) -> Vec<(String, NodeVec<DenseStorage<CopyDensity>>)> {
        let mut ret = Vec::new();

        for (em_id, m_step_result) in self.m.iter().enumerate() {
            for (m_id, m_step_once_result) in m_step_result.iter().enumerate() {
                match m_step_once_result {
                    MStepResult::Update(dbg, _) => {
                        let copy_densities = dbg.to_node_copy_densities();
                        let label = format!("{}#{}", em_id, m_id);
                        ret.push((label, copy_densities));
                    }
                    _ => {}
                }
            }
        }
        ret
    }
}

//
// For M-step
//

#[derive(Clone)]
pub enum MStepResult<K: KmerLike> {
    /// Initial dbg
    Init(FloatDbg<K>),
    /// The found negative cycle improved q-score
    Update(FloatDbg<K>, QScore),
    /// A negative cycle was found, but it did not improve q-score
    NoImprove(FloatDbg<K>, QScore),
    /// Any negative cycle was not found.
    NoNegCycle,
}
