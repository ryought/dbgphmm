use serde::ser::{SerializeSeq, Serializer};
use serde::Serialize;
use serde_json::ser::{Formatter, PrettyFormatter};
use std::io::{self, Write};

/// CustomFormatter (1): EscapeSlashes
struct EscapeSlashes;

impl Formatter for EscapeSlashes {
    fn write_string_fragment<W: ?Sized + Write>(
        &mut self,
        writer: &mut W,
        fragment: &str,
    ) -> io::Result<()> {
        writer.write_all(fragment.replace('/', "\\/").as_bytes())
    }
}

/// CustomFormatter (2): Rows
#[derive(Default)]
struct RowsFormatter {
    // pretty: PrettyFormatter<'static>,
    depth: usize,
}

impl RowsFormatter {
    pub fn new() -> Self {
        RowsFormatter { depth: 0 }
    }
}

impl Formatter for RowsFormatter {
    fn begin_array<W: ?Sized + Write>(&mut self, w: &mut W) -> io::Result<()> {
        self.depth += 1;
        w.write_all(b"[")
    }
    fn end_array<W: ?Sized + Write>(&mut self, w: &mut W) -> io::Result<()> {
        self.depth -= 1;
        if self.depth == 0 {
            w.write_all(b"\n");
        }
        w.write_all(b"]")
    }
    fn begin_array_value<W: ?Sized + Write>(&mut self, w: &mut W, first: bool) -> io::Result<()> {
        if !first {
            w.write_all(b",");
        }
        if self.depth == 1 {
            w.write_all(b"\n\t");
        }
        Ok(())
    }
    fn end_array_value<W: ?Sized + Write>(&mut self, w: &mut W) -> io::Result<()> {
        Ok(())
    }
}

#[derive(Serialize)]
struct Hoge {
    s: String,
    a: Vec<usize>,
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn json_custom_pretty() {
        // ref: https://github.com/serde-rs/json/issues/345#issuecomment-636215611
        let rows: Vec<Hoge> = vec![
            Hoge {
                s: "hoge/11".to_owned(),
                a: vec![1, 2, 3, 4],
            },
            Hoge {
                s: "fuga/12".to_owned(),
                a: vec![1, 2, 3, 4],
            },
        ];
        let out = std::io::stdout();
        // let mut ser = serde_json::Serializer::with_formatter(out, EscapeSlashes);
        let mut ser = serde_json::Serializer::with_formatter(out, RowsFormatter::default());
        // let mut ser = serde_json::Serializer::pretty(out);

        // (1) batch
        // rows.serialize(&mut ser).unwrap();

        // (2) stream
        let mut seq = ser.serialize_seq(None).unwrap();
        for row in rows.iter() {
            seq.serialize_element(&row).unwrap();
        }
        seq.end();
    }
}
