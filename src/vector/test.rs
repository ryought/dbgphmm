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
}
