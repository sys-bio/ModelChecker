//! A report as JSON, for the web page. Written by hand to keep the core
//! free of dependencies.

use crate::Report;
use crate::codes;

pub fn report(r: &Report) -> String {
    let mut s = String::from("{");
    field_opt_num(&mut s, "level", r.level);
    s.push(',');
    field_opt_num(&mut s, "version", r.version);
    s.push(',');
    field_opt_str(&mut s, "modelId", r.model_id.as_deref());
    s.push(',');
    field_opt_str(&mut s, "modelName", r.model_name.as_deref());
    s.push_str(",\"findings\":[");
    for (i, f) in r.findings.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push('{');
        field_str(&mut s, "code", f.code);
        s.push(',');
        field_str(
            &mut s,
            "title",
            codes::lookup(f.code).map(|c| c.title).unwrap_or(""),
        );
        s.push(',');
        field_str(&mut s, "severity", f.severity.as_str());
        s.push(',');
        field_str(&mut s, "kind", f.element.kind());
        s.push(',');
        field_str(&mut s, "id", &f.element.id());
        s.push(',');
        field_str(&mut s, "message", &f.message);
        s.push(',');
        field_opt_str(&mut s, "witness", f.witness.as_deref());
        s.push('}');
    }
    s.push_str("]}");
    s
}

/// Every code with its title and explanation, as JSON.
pub fn codes() -> String {
    let mut s = String::from("[");
    for (i, c) in codes::ALL.iter().enumerate() {
        if i > 0 {
            s.push(',');
        }
        s.push('{');
        field_str(&mut s, "code", c.code);
        s.push(',');
        field_str(&mut s, "title", c.title);
        s.push(',');
        field_str(&mut s, "help", &squash(c.help));
        s.push('}');
    }
    s.push(']');
    s
}

/// Collapses the runs of spaces left by source line continuations.
fn squash(t: &str) -> String {
    t.split_whitespace().collect::<Vec<_>>().join(" ")
}

fn field_str(s: &mut String, k: &str, v: &str) {
    string(s, k);
    s.push(':');
    string(s, v);
}

fn field_opt_str(s: &mut String, k: &str, v: Option<&str>) {
    string(s, k);
    s.push(':');
    match v {
        Some(v) => string(s, v),
        None => s.push_str("null"),
    }
}

fn field_opt_num(s: &mut String, k: &str, v: Option<u32>) {
    string(s, k);
    s.push(':');
    match v {
        Some(v) => s.push_str(&v.to_string()),
        None => s.push_str("null"),
    }
}

fn string(s: &mut String, v: &str) {
    s.push('"');
    for c in v.chars() {
        match c {
            '"' => s.push_str("\\\""),
            '\\' => s.push_str("\\\\"),
            '\n' => s.push_str("\\n"),
            '\r' => s.push_str("\\r"),
            '\t' => s.push_str("\\t"),
            c if (c as u32) < 0x20 => s.push_str(&format!("\\u{:04x}", c as u32)),
            c => s.push(c),
        }
    }
    s.push('"');
}
