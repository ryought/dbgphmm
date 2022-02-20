//!
//! Trait for types which have a unit of addition/multiplication
//!
use crate::prob::Prob;

///
/// get a unit of addition
///
pub trait UnitAdd: Copy + PartialEq {
    fn unit_add() -> Self;
    fn is_unit_add(&self) -> bool {
        *self == Self::unit_add()
    }
}

impl UnitAdd for u32 {
    fn unit_add() -> u32 {
        0
    }
}
impl UnitAdd for usize {
    fn unit_add() -> usize {
        0
    }
}
impl UnitAdd for f64 {
    fn unit_add() -> f64 {
        0.0
    }
}
impl UnitAdd for Prob {
    fn unit_add() -> Prob {
        Prob::from_log_prob(f64::NEG_INFINITY)
    }
}

///
/// get a unit of multiplication
///
pub trait UnitMul: Copy + PartialEq {
    fn unit_mul() -> Self;
    fn is_unit_mul(&self) -> bool {
        *self == Self::unit_mul()
    }
}

impl UnitMul for u32 {
    fn unit_mul() -> u32 {
        1
    }
}
impl UnitMul for usize {
    fn unit_mul() -> usize {
        1
    }
}
impl UnitMul for f64 {
    fn unit_mul() -> f64 {
        1.0
    }
}
impl UnitMul for Prob {
    fn unit_mul() -> Prob {
        Prob::from_log_prob(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_u32() {
        let e = u32::unit_add();
        println!("{} {} {}", e, e.is_unit_add(), 10u32.is_unit_add());
    }
}
