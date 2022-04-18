use super::Emission;
use crate::common::VALID_BASES;
use crate::hmmv2::params::PHMMParams;
use crate::prob::Prob;
use rand::prelude::*;

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
    let base = *VALID_BASES
        .choose_weighted(rng, |_base| param.p_random.to_value())
        .unwrap();
    Emission::Base(base)
}

///
/// Pick an emission at `Match` state (with a state emission), according to the PHMM params.
/// It selects a state emission with probability `p_match`, and other emissions with probability
/// `p_mismatch`.
///
pub fn pick_match_emission<R: Rng>(rng: &mut R, emission: u8, param: &PHMMParams) -> Emission {
    let base = *VALID_BASES
        .choose_weighted(rng, |&base| {
            if base == emission {
                param.p_match.to_value()
            } else {
                param.p_mismatch.to_value()
            }
        })
        .unwrap();
    Emission::Base(base)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_xoshiro::Xoshiro256PlusPlus;

    #[test]
    fn picker_pick_with_prob() {
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(0);

        // (1) assert that does not pick p=0 element
        for _ in 0..10 {
            let picked = pick_with_prob(
                &mut rng,
                &[(b'a', Prob::from_prob(0.0)), (b'b', Prob::from_prob(1.0))],
            );
            assert_eq!(picked, b'b');
        }

        // (2) assert that p=0.5 two elements evenly
        for _ in 0..100 {
            let picked = pick_with_prob(
                &mut rng,
                &[(b'a', Prob::from_prob(0.5)), (b'b', Prob::from_prob(0.5))],
            );
            // println!("{:?}", picked);
        }
    }
}
