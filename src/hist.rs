//!
//! Histogram counter
//!
use fnv::FnvHashMap as HashMap;

///
/// Histogram counter struct
///
/// * new
/// * add
/// * get
/// * range
/// * iter
///
#[derive(Clone, Debug)]
pub struct Hist(HashMap<usize, usize>);

impl Hist {
    ///
    /// Create an empty histogram counter
    ///
    pub fn new() -> Self {
        Hist(HashMap::default())
    }
    ///
    /// Increment a count of the value
    ///
    pub fn add(&mut self, value: usize) {
        *self.0.entry(value).or_insert(0) += 1;
    }
    ///
    /// Get count of the value
    ///
    pub fn get(&self, value: usize) -> usize {
        match self.0.get(&value) {
            Some(&count) => count,
            None => 0,
        }
    }
    ///
    /// Get the min/max of the values
    ///
    pub fn range(&self) -> Option<(usize, usize)> {
        if !self.0.is_empty() {
            let min = *self.0.keys().min().unwrap();
            let max = *self.0.keys().max().unwrap();
            Some((min, max))
        } else {
            None
        }
    }
    ///
    /// Iterate over values and its counts
    ///
    pub fn iter(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.0.iter().map(|(k, v)| (*k, *v))
    }
    ///
    /// Get the number of elements stored in the histogram
    ///
    pub fn len(&self) -> usize {
        self.0.values().sum()
    }
}

impl std::fmt::Display for Hist {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self.range() {
            Some((min, max)) => (min..=max).try_for_each(|value| {
                write!(
                    f,
                    "{}:{}{}",
                    value,
                    self.get(value),
                    if value != max { "," } else { "" }
                )
            }),
            None => Ok(()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn hist_test() {
        let mut h = Hist::new();
        h.add(10);
        h.add(10);
        h.add(2);
        h.add(10);
        h.add(2);
        h.add(3);
        println!("{:?}", h);
        assert_eq!(h.range(), Some((2, 10)));
        assert_eq!(h.get(10), 3);
        assert_eq!(h.get(2), 2);
        assert_eq!(h.get(3), 1);
        assert_eq!(h.get(0), 0);
        let c: Vec<_> = h.iter().collect();
        assert_eq!(c, vec![(2, 2), (3, 1), (10, 3)]);
        println!("{}", h);
        assert_eq!(h.to_string(), "2:2,3:1,4:0,5:0,6:0,7:0,8:0,9:0,10:3");
    }
}
