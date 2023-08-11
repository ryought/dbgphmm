//!
//! f64 that implements Ord to sort
//!

#[derive(Clone, Debug, Copy, PartialOrd, PartialEq)]
pub struct NonNanF64(f64);

impl NonNanF64 {
    pub fn new(value: f64) -> NonNanF64 {
        if value.is_nan() {
            panic!("NonNanF64 has Nan");
        }
        NonNanF64(value)
    }
}

impl From<f64> for NonNanF64 {
    fn from(value: f64) -> Self {
        NonNanF64::new(value)
    }
}

impl Eq for NonNanF64 {}

impl Ord for NonNanF64 {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.0.partial_cmp(&other.0).unwrap()
    }
}
