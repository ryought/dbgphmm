/// https://doc.rust-lang.org/rust-by-example/generics/new_types.html
/// new type idiom
struct Years(i64);
struct Days(i64);
impl Years {
    pub fn to_days(&self) -> Days {
        Days(self.0 * 365)
    }
}
impl Days {
    pub fn to_years(&self) -> Years {
        Years(self.0 / 365)
    }
}
fn summary_years(years: &Years) -> String {
    if years.0 > 10 {
        return "hoge".to_string();
    } else {
        return "fuga".to_string();
    }
}

pub fn test() {
    let y = Years(50);
    let d = y.to_days();
    println!("{}", summary_years(&y));
}
