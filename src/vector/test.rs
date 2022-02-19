#[cfg(test)]
mod tests {
    use super::super::{DenseStorage, SparseStorage, Storage, Vector};
    use petgraph::graph::NodeIndex;

    #[test]
    fn dense_vector() {
        let mut v: Vector<DenseStorage<u32>> = Vector::new(5, 0);
        v[0] = 100;
        v[3] = 222;
        assert_eq!(v[0], 100);
        assert_eq!(v[1], 0);
        let w: Vec<(usize, u32)> = v.iter().collect();
        println!("{:?}", v);
        println!("{:?}", w);
        assert_eq!(w, vec![(0, 100), (1, 0), (2, 0), (3, 222), (4, 0),]);
    }

    #[test]
    fn dense_vector_add_mul() {
        let mut v1: Vector<DenseStorage<u32>> = Vector::new(4, 0);
        v1[0] = 120;
        v1[3] = 111;
        let mut v2: Vector<DenseStorage<u32>> = Vector::new(4, 0);
        v2[0] = 1;
        v2[2] = 111;
        v2[3] = 1;
        // v1 + v2
        let added = &v1 + &v2;
        let muled = &v1 * &v2;
        println!("{:?}", added);
        assert_eq!(added[0], 120 + 1);
        assert_eq!(added[1], 0 + 0);
        assert_eq!(added[2], 0 + 111);
        assert_eq!(added[3], 111 + 1);
        println!("{:?}", muled);
        assert_eq!(muled[0], 120 * 1);
        assert_eq!(muled[1], 0 * 0);
        assert_eq!(muled[2], 0 * 111);
        assert_eq!(muled[3], 111 * 1);
        // v1 + v2
        v1 += &v2;
        assert_eq!(v1[0], 120 + 1);
        assert_eq!(v1[1], 0 + 0);
        assert_eq!(v1[2], 0 + 111);
        assert_eq!(v1[3], 111 + 1);
    }

    #[test]
    fn dense_vector_add_mul_constant() {
        let mut v1: Vector<DenseStorage<u32>> = Vector::new(4, 0);
        v1[0] = 120;
        v1[3] = 111;
        let added = v1 + 10;
        assert_eq!(added[0], 120 + 10);
        assert_eq!(added[1], 10);
        assert_eq!(added[2], 10);
        assert_eq!(added[3], 111 + 10);

        let mut v2: Vector<DenseStorage<u32>> = Vector::new(4, 0);
        v2[0] = 16;
        v2[2] = 77;
        let muled = v2 * 10;
        assert_eq!(muled[0], 160);
        assert_eq!(muled[1], 0);
        assert_eq!(muled[2], 770);
        assert_eq!(muled[3], 0);
    }

    #[test]
    fn sparse_vector() {
        let mut v: Vector<SparseStorage<u32>> = Vector::new(5, 0);
        v[0] = 100;
        v[3] = 222;
        assert_eq!(v[0], 100);
        assert_eq!(v[1], 0);
        let w: Vec<(usize, u32)> = v.iter().collect();
        println!("{:?}", v);
        println!("{:?}", w);
        assert_eq!(w, vec![(0, 100), (3, 222)]);
    }

    #[test]
    fn node_vector() {
        let mut v: Vector<DenseStorage<u32>, NodeIndex> = Vector::new(5, 0);
        v[NodeIndex::new(0)] = 111;
        v[NodeIndex::new(3)] = 222;
        assert_eq!(v[NodeIndex::new(0)], 111);
        assert_eq!(v[NodeIndex::new(1)], 0);
        let w: Vec<(NodeIndex, u32)> = v.iter().collect();
        println!("{:?}", v);
        println!("{:?}", w);
        assert_eq!(
            w,
            vec![
                (NodeIndex::new(0), 111),
                (NodeIndex::new(1), 0),
                (NodeIndex::new(2), 0),
                (NodeIndex::new(3), 222),
                (NodeIndex::new(4), 0),
            ]
        );
    }

    #[test]
    fn vector_conversion() {
        let default_value = 2;
        let mut v: Vector<DenseStorage<u32>, usize> = Vector::new(5, default_value);
        v[0] = 100;
        v[3] = 222;
        let w = v.to_sparse(default_value);
        println!("{:?}", v);
        println!("{:?}", w);
        assert_eq!(v[0], w[0]);
        assert_eq!(v[1], w[1]);
        assert_eq!(v[2], w[2]);
        assert_eq!(v[3], w[3]);
        assert_eq!(v[4], w[4]);
        let v2 = v.to_dense();
        println!("{:?}", v2);
        assert_eq!(v[0], v2[0]);
        assert_eq!(v[1], v2[1]);
        assert_eq!(v[2], v2[2]);
        assert_eq!(v[3], v2[3]);
        assert_eq!(v[4], v2[4]);
    }

    #[test]
    fn vector_ops_between_different_storage_type() {
        let mut v1: Vector<DenseStorage<u32>, usize> = Vector::new(5, 0);
        v1[0] = 33;
        v1[3] = 55;
        let mut v2: Vector<SparseStorage<u32>, usize> = Vector::new(5, 0);
        v2[0] = 22;
        v2[2] = 11;
        println!("{:?}", v1);
        println!("{:?}", v2);
        let w = &v1 + &v2;
        println!("{:?}", w);
        assert!(w.is_dense());
        assert_eq!(w[0], 33 + 22);
        assert_eq!(w[1], 0);
        assert_eq!(w[2], 11);
        assert_eq!(w[3], 55);
    }
    #[test]
    fn vector_to_sparse_by_indexes() {
        let mut v: Vector<DenseStorage<u32>, usize> = Vector::new(10, 0);
        v[0] = 33;
        v[3] = 55;
        v[5] = 110;
        println!("{:?}", v);
        let w = v.to_sparse_by_indexes(0, &[5, 3]);
        println!("{:?}", w);
        let e: Vec<(usize, u32)> = w.iter().collect();
        assert_eq!(e, vec![(5, 110), (3, 55)]);
    }
    #[test]
    fn vector_approx_eq() {
        // dense
        let mut v: Vector<DenseStorage<f64>, usize> = Vector::new(10, 0.0);
        v[0] = 99.9;
        v[3] = 82.2;
        let mut w: Vector<DenseStorage<f64>, usize> = Vector::new(10, 0.0);
        w[0] = 99.899999;
        w[3] = 82.2;
        assert!(!abs_diff_eq!(v, w));
        assert!(abs_diff_eq!(v, w, epsilon = 0.1));
        let mut w2: Vector<DenseStorage<f64>, usize> = Vector::new(10, 0.0);
        w2[0] = 99.9;
        w2[3] = 82.2;
        assert!(abs_diff_eq!(v, w2));

        // sparse
        let mut v: Vector<SparseStorage<f64>, usize> = Vector::new(10, 0.0);
        v[0] = 99.9;
        v[3] = 82.2;
        let mut w: Vector<SparseStorage<f64>, usize> = Vector::new(10, 0.0);
        w[0] = 99.899999;
        w[3] = 82.2;
        assert!(!abs_diff_eq!(v, w));
        assert!(abs_diff_eq!(v, w, epsilon = 0.1));
    }
    #[test]
    fn vector_add_mul_more() {
        // add
        // Dense([10, 3, 3, 3, 10])
        let mut a: Vector<DenseStorage<u32>, usize> = Vector::from_vec(5, 3, &[(0, 10), (4, 10)]);
        // Sparse([2, 1, 1, 1, 1])
        let mut b: Vector<SparseStorage<u32>, usize> = Vector::from_vec(5, 1, &[(0, 2)]);
        assert_eq!((&a + &a).to_vec(), vec![20, 6, 6, 6, 20]);
        assert_eq!((&a + &b).to_vec(), vec![12, 4, 4, 4, 11]);
        assert_eq!((&b + &a).to_vec(), vec![12, 4, 4, 4, 11]);
        assert_eq!((&b + &b).to_vec(), vec![4, 2, 2, 2, 2]);

        // mul
        // Dense([10, 0, 1, 0, 5])
        let mut a: Vector<DenseStorage<u32>, usize> =
            Vector::from_vec(5, 3, &[(0, 10), (1, 0), (2, 1), (3, 0), (4, 5)]);
        // Sparse([0, 1, 10, 0, 1])
        let mut b: Vector<SparseStorage<u32>, usize> =
            Vector::from_vec(5, 0, &[(1, 1), (2, 10), (4, 1)]);
        assert_eq!((&a * &a).to_vec(), vec![100, 0, 1, 0, 25]);
        assert_eq!((&a * &b).to_vec(), vec![0, 0, 10, 0, 5]);
        assert_eq!((&b * &a).to_vec(), vec![0, 0, 10, 0, 5]);
        assert_eq!((&b * &b).to_vec(), vec![0, 1, 100, 0, 1]);

        // mul
        // Dense([10, 0, 1, 2, 5])
        let mut a: Vector<DenseStorage<u32>, usize> =
            Vector::from_vec(5, 3, &[(0, 10), (1, 0), (2, 1), (3, 2), (4, 5)]);
        // Sparse([1, 1, 3, 1, 0])
        let mut b: Vector<SparseStorage<u32>, usize> = Vector::from_vec(5, 1, &[(2, 3), (4, 0)]);
        assert_eq!((&a * &a).to_vec(), vec![100, 0, 1, 4, 25]);
        assert_eq!((&a * &b).to_vec(), vec![10, 0, 3, 2, 0]);
        assert_eq!((&b * &a).to_vec(), vec![10, 0, 3, 2, 0]);
        assert_eq!((&b * &b).to_vec(), vec![1, 1, 9, 1, 0]);
    }
}
