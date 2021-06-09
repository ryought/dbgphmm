pub fn test2() {
    let mut xs: Vec<Vec<u8>> = Vec::new();
    for i in 1..10 {
        let mut x = Vec::new();
        for j in 1..5 {
            x.push(j + 2);
        }
        xs.push(x);
    }
    println!("xs: {:?}", xs);
}

/*
pub fn test() {
    let mut xs: Vec<&[u8]> = Vec::new();
    for i in 1..10 {
        let mut x = Vec::new();
        for j in 1..5 {
            x.push(j + 2);
        }
        xs.push(&x);
    }
    println!("xs: {:?}", xs);
}
*/
