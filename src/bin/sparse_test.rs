use dbgphmm::veclike::{SparseVec, VecLike};
fn main() {
    let mut v: SparseVec<usize> = SparseVec::new(10, 0);
    v.set(5, 101);
    v.set(0, 1111);
    v.set(3, 11);
    v.set(7, 89);
    for x in v.iter_sparse() {
        println!("{:?}", x);
    }
}
