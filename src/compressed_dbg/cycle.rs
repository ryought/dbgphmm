//!
//! Cycle-related implementations of CompressedDBG
//!
use super::CompressedDBG;
use crate::cycles::CycleDirection;
use crate::graph::Node;
use itertools::iproduct;

impl CompressedDBG {
    pub fn n_cycles(&self) -> usize {
        self.cycles.len()
    }
    pub fn cycle_components(&self, cycle_id: usize) -> &[(Node, CycleDirection)] {
        &self.cycles[cycle_id]
    }
    /// cycle is a sequence of (node, direction)
    /// when is_up=true, node with direction=Forward will get +1.
    /// all nodes should be >=0 after this update
    pub fn is_acceptable(&self, copy_nums: &[u32], cycle_id: usize, is_up: bool) -> bool {
        self.cycle_components(cycle_id).iter().all(|(v, dir)| {
            if is_up && *dir == CycleDirection::Reverse {
                copy_nums[v.0] > 0
            } else if !is_up && *dir == CycleDirection::Forward {
                copy_nums[v.0] > 0
            } else {
                true
            }
        })
    }
    pub fn update_by_cycle(&self, copy_nums: &[u32], cycle_id: usize, is_up: bool) -> Vec<u32> {
        let mut new_copy_nums = copy_nums.to_vec();
        for (v, dir) in self.cycle_components(cycle_id).iter() {
            if is_up {
                match dir {
                    CycleDirection::Forward => new_copy_nums[v.0] += 1,
                    CycleDirection::Reverse => new_copy_nums[v.0] -= 1,
                }
            } else {
                match dir {
                    CycleDirection::Forward => new_copy_nums[v.0] -= 1,
                    CycleDirection::Reverse => new_copy_nums[v.0] += 1,
                }
            }
        }
        new_copy_nums
    }
    pub fn update_cycle_vec_by_cycle(
        &self,
        cycle_vec: &[u32],
        cycle_id: usize,
        is_up: bool,
    ) -> Vec<u32> {
        let mut new_cycle_vec = cycle_vec.to_vec();
        if is_up {
            new_cycle_vec[cycle_id] += 1;
        } else {
            new_cycle_vec[cycle_id] -= 1;
        }
        new_cycle_vec
    }
    /// Convert copy_nums into cycle_vec
    /// For each cycle basis,
    /// (the amount of cycle basis in the cycle vec) == (copy_num of cycle key edge)
    /// and the total amount will be the cycle vec
    pub fn cycle_vec_from_copy_nums(&self, copy_nums: &[u32]) -> Vec<u32> {
        let mut cycle_vec: Vec<i32> = vec![0; self.n_cycles()];
        for cycle_id in 0..self.n_cycles() {
            let (node, _) = self.cycle_components(cycle_id).first().unwrap();
            let count = copy_nums[node.0];
            cycle_vec[cycle_id] += count as i32;
        }
        cycle_vec
            .iter()
            .map(|&x| {
                assert!(x >= 0);
                x as u32
            })
            .collect()
    }
    /// Convert cycle_vec into copy_nums
    /// Inverse function of cycle_vec_from_copy_nums
    pub fn copy_nums_from_cycle_vec(&self, cycle_vec: &[u32]) -> Vec<u32> {
        let mut copy_nums: Vec<i32> = vec![0; self.n_kmers()];
        for cycle_id in 0..self.n_cycles() {
            let count = cycle_vec[cycle_id];
            for (node, dir) in self.cycle_components(cycle_id).iter() {
                match dir {
                    CycleDirection::Forward => copy_nums[node.0] += count as i32,
                    CycleDirection::Reverse => copy_nums[node.0] -= count as i32,
                }
            }
        }
        copy_nums
            .iter()
            .map(|&x| {
                assert!(x >= 0);
                x as u32
            })
            .collect()
    }
    pub fn is_all_cycle_consistent(&self, init_copy_nums: &[u32]) {
        for i in 0..self.n_cycles() {
            let x = self.is_acceptable(init_copy_nums, i, true);
            let y = self.is_acceptable(init_copy_nums, i, false);
            println!("i={}, x={}, y={}", i, x, y);
        }
    }
    pub fn validate_cycles(&self, copy_nums_true: &[u32]) {
        for i in 0..self.n_cycles() {
            let is_a = self.is_acceptable(&copy_nums_true, i, true);
            let is_b = self.is_acceptable(&copy_nums_true, i, false);
            assert!(is_a);
            assert!(is_b);
            let new_a = self.update_by_cycle(&copy_nums_true, i, true);
            let new_b = self.update_by_cycle(&copy_nums_true, i, false);
            assert!(self.is_consistent_copy_num(&new_a));
            assert!(self.is_consistent_copy_num(&new_b));
        }
    }
    pub fn cycle_and_direction_candidates(&self, copy_nums: &[u32]) -> Vec<(usize, bool)> {
        let candidates: Vec<(usize, bool)> = iproduct!((0..self.n_cycles()), &[true, false])
            .filter(|(cycle_id, &is_up)| self.is_acceptable(copy_nums, *cycle_id, is_up))
            .map(|(cycle_id, is_up)| (cycle_id, *is_up))
            .collect();
        candidates
    }
}
