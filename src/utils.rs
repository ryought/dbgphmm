pub fn is_equal_as_set<T: PartialEq>(xs: &[T], ys: &[T]) -> bool {
    xs.iter().all(|x| ys.contains(x)) && ys.iter().all(|y| xs.contains(y))
}

///
/// Check if all items from iterator is the same value
/// If so, return the unique value.
/// Otherwise (including the iterator has no value) return none.
///
pub fn all_same_value<T: PartialEq + Clone, I: Iterator<Item = T>>(mut iter: I) -> Option<T> {
    match iter.next() {
        Some(first) => {
            if iter.all(|item| item == first) {
                Some(first)
            } else {
                None
            }
        }
        None => None,
    }
}

///
/// TODO
/// * change return type from Vec<T> into Option<Vec<T>>?
///
pub fn unwrap_all<T>(options: Vec<Option<T>>) -> Vec<T> {
    options
        .into_iter()
        .map(|option| option.expect(""))
        .collect()
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

    #[test]
    fn all_same_value_test() {
        {
            let xs = vec![0, 1, 2, 5, 3];
            assert_eq!(all_same_value(xs.iter()), None);
        }
        {
            let xs = vec![0, 0, 0];
            assert_eq!(all_same_value(xs.iter()), Some(&0));
        }
        {
            let xs = vec![0, 0, 1];
            assert_eq!(all_same_value(xs.iter()), None);
        }
        {
            let xs = vec![5];
            assert_eq!(all_same_value(xs.iter()), Some(&5));
        }
        {
            let xs: Vec<usize> = Vec::new();
            assert_eq!(all_same_value(xs.iter()), None);
        }
    }
}
