use std::fs::File;
use std::io::Write;

fn f<F: Write>(f: &mut F) {
    writeln!(f, "hoge");
}

fn main() {
    let is_stdout = true;
    if is_stdout {
        let mut file = File::create("hoge.txt").unwrap();
        f(&mut file);
    } else {
        f(&mut std::io::stdout());
    }
}
