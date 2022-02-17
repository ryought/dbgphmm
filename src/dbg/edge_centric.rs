//!
//! Edge-centric de Bruijn graph implementations
//!
//! * Node = k-1-mer overlap
//! * Edge = k-mer
//!
//! # Features
//!
//! * visualization by cytoscape
//! * min-flow optimization
//! * quick-check of copy number consistency
//!
use crate::common::CopyNum;
use crate::kmer::kmer::{Kmer, KmerLike};
use petgraph::graph::{DiGraph, EdgeIndex, NodeIndex};
pub mod impls;

///
/// (Edge-centric) De bruijn graph struct
///
pub struct EDbg<N: EDbgNode, E: EDbgEdge> {
    k: usize,
    pub graph: DiGraph<N, E>,
}

///
/// Trait for nodes in edge-centric dbg `EDbg`
///
pub trait EDbgNode {
    type Kmer: KmerLike;
    ///
    /// k-1-mer of this node of the EDbg
    fn km1mer(&self) -> &Self::Kmer;
    ///
    /// Copy number count of this node in EDbg
    fn copy_num(&self) -> CopyNum;
}

///
/// Trait for edges in edge-centric dbg `EDbg`
///
pub trait EDbgEdge {
    type Kmer: KmerLike;
    ///
    /// k-mer of this edge of the EDbg
    fn kmer(&self) -> &Self::Kmer;
    ///
    /// Copy number count of this edge in EDbg
    fn copy_num(&self) -> CopyNum;
}
