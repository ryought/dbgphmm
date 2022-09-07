//!
//! Q score calculation functions
//! after forward/backward calculation
//!
//! - q_score
//! - q_score_diff
//!

#[derive(Clone, Debug, Copy, Default)]
pub struct QScore {
    /// init score
    pub init: f64,
    /// trans score
    pub trans: f64,
    /// prior score
    pub prior: f64,
}

impl QScore {
    ///
    /// Constructor
    ///
    pub fn new(init: f64, trans: f64, prior: f64) -> Self {
        QScore { init, trans, prior }
    }
    ///
    /// total score
    ///
    pub fn total(&self) -> f64 {
        self.init + self.trans + self.prior
    }
}

impl std::fmt::Display for QScore {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{}(init={} trans={} prior={})",
            self.total(),
            self.init,
            self.trans,
            self.prior
        )
    }
}
