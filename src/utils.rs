pub fn is_equal_as_set<T: PartialEq>(xs: &[T], ys: &[T]) -> bool {
    xs.iter().all(|x| ys.contains(x)) && ys.iter().all(|y| xs.contains(y))
}

//
// tests
//

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_equal_as_set_test_1() {
        let xs = vec![0, 1, 2, 5, 3];
        let ys = vec![1, 5, 0, 3, 2];
        assert!(is_equal_as_set(&xs, &ys));

        let zs = vec![1, 5, 3, 2];
        assert!(!is_equal_as_set(&xs, &zs));
    }
}
