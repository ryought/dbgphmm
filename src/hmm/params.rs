use crate::prob::Prob;

#[derive(Debug, Clone, PartialEq)]
pub struct PHMMParams {
    pub p_mismatch: Prob,
    pub p_match: Prob,
    pub p_random: Prob,
    pub p_gap_open: Prob,
    pub p_gap_ext: Prob,
    pub p_end: Prob,
    pub p_MM: Prob,
    pub p_IM: Prob,
    pub p_DM: Prob,
    pub p_MI: Prob,
    pub p_II: Prob,
    pub p_DI: Prob,
    pub p_MD: Prob,
    pub p_ID: Prob,
    pub p_DD: Prob,
    pub n_max_gaps: u32,
}

impl PHMMParams {
    pub fn new(
        p_mismatch: Prob,
        p_gap_open: Prob,
        p_gap_ext: Prob,
        p_end: Prob,
        n_max_gaps: u32,
    ) -> PHMMParams {
        PHMMParams {
            p_mismatch,
            p_gap_open,
            p_gap_ext,
            p_end,
            p_DD: p_gap_ext,
            p_II: p_gap_ext,
            p_MI: p_gap_open,
            p_MD: p_gap_open,
            p_ID: p_gap_open,
            p_DI: p_gap_open,
            // p_MM: 1 - p_gap_open - p_gap_open,
            p_MM: Prob::from_prob(1.0 - 2.0 * p_gap_open.to_value() - p_end.to_value()),
            // p_DM: 1 - p_gap_open - p_gap_ext,
            p_DM: Prob::from_prob(
                1.0 - p_gap_open.to_value() - p_gap_ext.to_value() - p_end.to_value(),
            ),
            p_IM: Prob::from_prob(
                1.0 - p_gap_open.to_value() - p_gap_ext.to_value() - p_end.to_value(),
            ),
            // p_match: 1 - p_mismatch
            p_match: Prob::from_prob(1.0 - p_mismatch.to_value()),
            p_random: Prob::from_prob(0.25),
            n_max_gaps,
        }
    }
    pub fn default() -> PHMMParams {
        PHMMParams::new(
            Prob::from_prob(0.01),    // mismatch
            Prob::from_prob(0.01),    // gap_open
            Prob::from_prob(0.01),    // gap_ext
            Prob::from_prob(0.00001), // end
            3,
        )
    }
}

impl std::fmt::Display for PHMMParams {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        writeln!(f, "p_mismatch: {}", self.p_mismatch);
        writeln!(f, "p_match: {}", self.p_match);
        writeln!(f, "p_gap_open: {}", self.p_gap_open);
        writeln!(f, "p_gap_ext: {}", self.p_gap_ext);
        writeln!(f, "p_end: {}", self.p_end);
        writeln!(f, "p_MM: {}", self.p_MM);
        writeln!(f, "p_IM: {}", self.p_IM);
        writeln!(f, "p_DM: {}", self.p_DM);
        writeln!(f, "p_MI: {}", self.p_MI);
        writeln!(f, "p_II: {}", self.p_II);
        writeln!(f, "p_DI: {}", self.p_DI);
        writeln!(f, "p_MD: {}", self.p_MD);
        writeln!(f, "p_ID: {}", self.p_ID);
        writeln!(f, "p_DD: {}", self.p_DD);
        writeln!(f, "n_max_gaps: {}", self.n_max_gaps);
        Ok(())
    }
}
