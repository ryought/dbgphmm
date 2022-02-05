use super::Emission;
use crate::hmm::params::PHMMParams;
use crate::prob::Prob;
use rand::prelude::*;

///
/// Array of valid DNA bases
///
const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

///
/// pick randomly from the choices with its own probability.
///
pub fn pick_with_prob<R: Rng, T: Copy>(rng: &mut R, choices: &[(T, Prob)]) -> T {
    choices
        .choose_weighted(rng, |item| item.1.to_value())
        .unwrap()
        .0
}

///
/// Pick an emission at `Ins` state, according to the PHMM params.
/// It selects uniformaly random on `ACGT`
///
pub fn pick_ins_emission<R: Rng>(rng: &mut R, param: &PHMMParams) -> Emission {
    let base = *BASES
        .choose_weighted(rng, |_base| param.p_random.to_value())
        .unwrap();
    Emission(Some(base))
}

///
/// Pick an emission at `Match` state (with a state emission), according to the PHMM params.
/// It selects a state emission with probability `p_match`, and other emissions with probability
/// `p_mismatch`.
///
pub fn pick_match_emission<R: Rng>(rng: &mut R, emission: u8, param: &PHMMParams) -> Emission {
    let base = *BASES
        .choose_weighted(rng, |&base| {
            if base == emission {
                param.p_match.to_value()
            } else {
                param.p_mismatch.to_value()
            }
        })
        .unwrap();
    Emission(Some(base))
}
