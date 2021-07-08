use serde::{Deserialize, Serialize};
use serde_json;

#[derive(Serialize, Deserialize, Debug, Default)]
pub struct DbgStats {
    pub k: usize,
    pub n_kmers: u32,
    pub n_starting_kmers: u32,
}

#[derive(Serialize, Deserialize, Debug, Default)]
pub struct CopyNumStats {
    pub total: u32,
    pub total_emitable: u32,
    pub is_consistent: bool,
    pub average: f32,
    pub min: u32,
    pub max: u32,
    pub n_zero_copy_kmer: u32,
    pub n_nonzero_copy_kmer: u32,
}

#[derive(Serialize, Deserialize, Debug, Default)]
pub struct DegreeStats {
    pub in_degs: [u32; 6],
    pub out_degs: [u32; 6],
}

#[derive(Serialize, Deserialize, Debug, Default)]
/// TODO needs unionfind?
pub struct PathStats {
    pub n_simple_paths: u32,
}

#[derive(Serialize, Deserialize, Debug, Default)]
pub struct CycleSummaryStats {
    pub n_cycles: u32,
    pub min_cycle_len: u32,
    pub max_cycle_len: u32,
    pub average_cycle_len: f32,
}

#[derive(Serialize, Deserialize, Debug, Default)]
pub struct CycleStats {
    pub id: usize,
    pub len: usize,
    pub n_reverse: usize,
}

#[derive(Serialize, Deserialize, Debug, Default)]
pub struct CompareResult {
    pub n_shared: u32,
    pub n_self_only: u32,
    pub n_other_only: u32,
}

#[derive(Serialize, Deserialize, Debug, Default)]
pub struct AllStats {
    dbg_stats: Option<DbgStats>,
    copy_num_stats: Option<CopyNumStats>,
}

pub fn test() {
    // let point = DbgStats::default();
    let point = AllStats {
        dbg_stats: None,
        copy_num_stats: Some(CopyNumStats::default()),
    };
    let serialized = serde_json::to_string_pretty(&point).unwrap();
    println!("{}", serialized);
}
