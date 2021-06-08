/// phmm dbg core
/// depends dbg for
/// - iter on kmers
/// - transitions (childs(v))
/// - transition probabilities (p(v, w))
/// - emissions (emit(v))
/// - copy nums (init(v))
///
/// forward
/// - init() -> F[0]
/// - step_emit_states(F[i-1]) -> F[i]
/// - step_silent_states(F[i]) -> F[i]
struct PHMM {}

impl PHMM {
    fn forward() -> Vec<layer> {}
    fn backward() -> Vec<layer> {}
    fn step(layer) -> layer {}
}
