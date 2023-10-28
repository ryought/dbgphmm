//!
//! Profile HMM calculation
//!
//! # Overview of calculation
//!
//! x = x[0],...,x[n-1] : Emissions of length n
//!
//! Forward
//! F[i][node]
//!  = P(emits x[0:i]=x[0],...,x[i-1] and ends at node) for 0<=i<=n
//!
//! Backward
//! B[i][node]
//!  = P(emits x[i:n] | starts from node) for 0<=i<=n
//!
//! Freq / Mapping
//! S[i][node]
//!  = P(node emits x[i])
//!  = P(emits x[0:i+1]=x[0],...,x[i] and ends at node) P(starts from node and then emits x[i+1:n]=x[i+1],...,x[n-1])
//!  = F[i+1] B[i+1]
//!
//! S[i] = S.init_table (i=-1; S[-1]=F[0]B[0])
//!        S.tables[i]  (0<=i<n)
//!
//! F[i] = F.init_table  (i=0; no emission)
//!        F.tables[i-1] (0<i<=n)
//!
//! B[i] = B.init_table  (i=n; no emission)
//!        B.tables[i]   (0<=i<n)
//!
pub mod active_nodes;
pub mod backward;
pub mod common;
pub mod forward;
pub mod freq;
pub mod hint;
pub mod mocks;
pub mod new;
pub mod params;
pub mod q;
pub mod sample;
pub mod speed;
pub mod table;
pub mod tests;
pub mod trans_table;

pub use freq::NodeFreqs;
pub use trans_table::EdgeFreqs;
