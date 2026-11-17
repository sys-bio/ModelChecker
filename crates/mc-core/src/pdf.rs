//! A report as a PDF, for the web page's "Save report as PDF". Written by
//! hand, like the JSON, to keep the core free of dependencies: the PDF uses
//! the standard fonts every reader has (Helvetica, Helvetica-Bold, Courier),
//! so nothing is embedded, and text is encoded as WinAnsi.

use std::collections::HashMap;

use crate::codes;
use crate::{Finding, Report, Severity};

/// What the report says about where the model came from and when it was
/// checked; the page knows these, the checker does not.
pub struct Meta<'a> {
    /// The file the model was loaded from ("pasted.ant" when pasted).
    pub file_name: &'a str,
    /// What was checked: "SBML Level 3 Version 2", "Antimony", ...
    pub source: &'a str,
    /// When the report was made, as the reader should see it.
    pub date: &'a str,
    /// The checker and its version.
    pub checker: &'a str,
}

// US Letter, in points.
const PAGE_W: f32 = 612.0;
const PAGE_H: f32 = 792.0;
const LEFT: f32 = 54.0;
const RIGHT: f32 = PAGE_W - 54.0;
const TOP: f32 = PAGE_H - 54.0;
/// Lowest baseline for body text; the footer sits below it.
const BOTTOM: f32 = 66.0;

type Rgb = (f32, f32, f32);
const TEXT: Rgb = (0.10, 0.11, 0.13);
const MUTED: Rgb = (0.40, 0.42, 0.46);
const RULE: Rgb = (0.82, 0.83, 0.86);
const ERROR: Rgb = (0.75, 0.11, 0.11);
const WARNING: Rgb = (0.70, 0.42, 0.00);
const INFO: Rgb = (0.20, 0.38, 0.62);

fn colour(s: Severity) -> Rgb {
    match s {
        Severity::Error => ERROR,
        Severity::Warning => WARNING,
        Severity::Info => INFO,
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Font {
    Regular,
    Bold,
    Mono,
}

impl Font {
    fn name(self) -> &'static str {
        match self {
            Font::Regular => "F1",
            Font::Bold => "F2",
            Font::Mono => "F3",
        }
    }
}

/// Advance widths (1/1000 em) of characters 32 to 126, from the fonts'
/// metrics files.
const HELVETICA: [u16; 95] = [
    278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278,
    278, // space to /
    556, 556, 556, 556, 556, 556, 556, 556, 556, 556, // 0 to 9
    278, 278, 584, 584, 584, 556, 1015, // : to @
    667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833, // A to M
    722, 778, 667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, // N to Z
    278, 278, 278, 469, 556, 333, // [ to `
    556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500, 222, 833, // a to m
    556, 556, 556, 556, 333, 500, 278, 556, 500, 722, 500, 500, 500, // n to z
    334, 260, 334, 584, // { to ~
];
const HELVETICA_BOLD: [u16; 95] = [
    278, 333, 474, 556, 556, 889, 722, 238, 333, 333, 389, 584, 278, 333, 278,
    278, // space to /
    556, 556, 556, 556, 556, 556, 556, 556, 556, 556, // 0 to 9
    333, 333, 584, 584, 584, 611, 975, // : to @
    722, 722, 722, 722, 667, 611, 778, 722, 278, 556, 722, 611, 833, // A to M
    722, 778, 667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611, // N to Z
    333, 278, 333, 584, 556, 333, // [ to `
    556, 611, 556, 611, 556, 333, 611, 611, 278, 278, 556, 278, 889, // a to m
    611, 611, 611, 611, 389, 556, 333, 611, 556, 778, 556, 556, 500, // n to z
    389, 280, 389, 584, // { to ~
];

fn char_width(font: Font, b: u8) -> u16 {
    match font {
        Font::Mono => 600,
        Font::Regular if (32..=126).contains(&b) => HELVETICA[(b - 32) as usize],
        Font::Bold if (32..=126).contains(&b) => HELVETICA_BOLD[(b - 32) as usize],
        // Accented letters and symbols: about the width of a letter.
        _ => 556,
    }
}

fn width(font: Font, size: f32, s: &[u8]) -> f32 {
    s.iter().map(|&b| char_width(font, b) as f32).sum::<f32>() * size / 1000.0
}

/// Text as WinAnsi bytes. Characters WinAnsi lacks are written as their
/// nearest ASCII ("->" for an arrow), or "?". Newlines are kept.
fn encode(s: &str) -> Vec<u8> {
    let mut out = Vec::with_capacity(s.len());
    for c in s.chars() {
        let code = c as u32;
        match c {
            '\n' => out.push(b'\n'),
            '\t' => out.extend_from_slice(b"  "),
            ' '..='~' => out.push(code as u8),
            '\u{a0}'..='\u{ff}' => out.push(code as u8),
            _ => {
                if let Some(b) = win_ansi(c) {
                    out.push(b);
                } else {
                    out.extend_from_slice(ascii_for(c).as_bytes());
                }
            }
        }
    }
    out
}

/// WinAnsi's characters in 0x80 to 0x9F.
fn win_ansi(c: char) -> Option<u8> {
    Some(match c {
        '\u{20ac}' => 0x80,
        '\u{201a}' => 0x82,
        '\u{0192}' => 0x83,
        '\u{201e}' => 0x84,
        '\u{2026}' => 0x85,
        '\u{2020}' => 0x86,
        '\u{2021}' => 0x87,
        '\u{02c6}' => 0x88,
        '\u{2030}' => 0x89,
        '\u{0160}' => 0x8a,
        '\u{2039}' => 0x8b,
        '\u{0152}' => 0x8c,
        '\u{017d}' => 0x8e,
        '\u{2018}' => 0x91,
        '\u{2019}' => 0x92,
        '\u{201c}' => 0x93,
        '\u{201d}' => 0x94,
        '\u{2022}' => 0x95,
        '\u{2013}' => 0x96,
        '\u{2014}' => 0x97,
        '\u{02dc}' => 0x98,
        '\u{2122}' => 0x99,
        '\u{0161}' => 0x9a,
        '\u{203a}' => 0x9b,
        '\u{0153}' => 0x9c,
        '\u{017e}' => 0x9e,
        '\u{0178}' => 0x9f,
        _ => return None,
    })
}

fn ascii_for(c: char) -> &'static str {
    match c {
        '\u{2192}' | '\u{27f6}' => "->",
        '\u{2190}' => "<-",
        '\u{21cc}' | '\u{2194}' => "<->",
        '\u{2264}' => "<=",
        '\u{2265}' => ">=",
        '\u{2260}' => "!=",
        '\u{2212}' => "-",
        '\u{2248}' => "~",
        '\u{221e}' => "inf",
        '\u{22c5}' | '\u{2219}' => "*",
        // Thin and other narrow spaces (dates and times use them).
        '\u{2002}'..='\u{200a}' | '\u{202f}' | '\u{205f}' => " ",
        _ => "?",
    }
}

/// Breaks text into lines at most `max` points wide, at spaces where it
/// can, else inside a word. Newlines start new lines.
fn wrap(font: Font, size: f32, max: f32, text: &[u8]) -> Vec<Vec<u8>> {
    let mut lines = Vec::new();
    for para in text.split(|&b| b == b'\n') {
        let mut line: Vec<u8> = Vec::new();
        for word in para.split(|&b| b == b' ') {
            let candidate = if line.is_empty() {
                word.to_vec()
            } else {
                [&line[..], b" ", word].concat()
            };
            if width(font, size, &candidate) <= max {
                line = candidate;
                continue;
            }
            if !line.is_empty() {
                lines.push(std::mem::take(&mut line));
            }
            // A word wider than the line is broken where it must be.
            let mut rest = word;
            while width(font, size, rest) > max {
                let mut n = 1;
                while n < rest.len() && width(font, size, &rest[..n + 1]) <= max {
                    n += 1;
                }
                lines.push(rest[..n].to_vec());
                rest = &rest[n..];
            }
            line = rest.to_vec();
        }
        lines.push(line);
    }
    lines
}

fn hex(s: &[u8]) -> String {
    let mut out = String::with_capacity(s.len() * 2 + 2);
    out.push('<');
    for b in s {
        out.push_str(&format!("{b:02X}"));
    }
    out.push('>');
    out
}

/// Pages being laid out, top to bottom.
struct Doc {
    pages: Vec<String>,
    page: String,
    /// The current baseline.
    y: f32,
}

impl Doc {
    fn new() -> Self {
        Doc {
            pages: Vec::new(),
            page: String::new(),
            y: TOP,
        }
    }

    fn new_page(&mut self) {
        self.pages.push(std::mem::take(&mut self.page));
        self.y = TOP;
    }

    /// Starts a new page unless this one is still empty.
    fn fresh_page(&mut self) {
        if !self.page.is_empty() {
            self.new_page();
        }
    }

    /// Starts a new page unless `h` points fit on this one.
    fn need(&mut self, h: f32) {
        if self.y - h < BOTTOM {
            self.new_page();
        }
    }

    fn text_at(&mut self, x: f32, y: f32, font: Font, size: f32, c: Rgb, s: &[u8]) {
        self.page.push_str(&format!(
            "BT /{} {size:.1} Tf {:.3} {:.3} {:.3} rg {x:.2} {y:.2} Td {} Tj ET\n",
            font.name(),
            c.0,
            c.1,
            c.2,
            hex(s)
        ));
    }

    fn rule(&mut self, y: f32) {
        self.page.push_str(&format!(
            "{:.3} {:.3} {:.3} RG 0.6 w {LEFT:.2} {y:.2} m {RIGHT:.2} {y:.2} l S\n",
            RULE.0, RULE.1, RULE.2
        ));
    }

    /// Moves down one line of `leading` and returns its baseline.
    fn next_line(&mut self, leading: f32) -> f32 {
        self.need(leading);
        self.y -= leading;
        self.y
    }

    /// Writes wrapped text from `x` to the right margin.
    #[allow(clippy::too_many_arguments)]
    fn para(&mut self, x: f32, font: Font, size: f32, c: Rgb, leading: f32, s: &str) {
        for line in wrap(font, size, RIGHT - x, &encode(s)) {
            let y = self.next_line(leading);
            self.text_at(x, y, font, size, c, &line);
        }
    }

    /// A label in grey followed by wrapped text, the label hanging.
    #[allow(clippy::too_many_arguments)]
    fn labelled(&mut self, label: &str, indent: f32, font: Font, size: f32, leading: f32, s: &str) {
        let lines = wrap(font, size, RIGHT - LEFT - indent, &encode(s));
        for (i, line) in lines.iter().enumerate() {
            let y = self.next_line(leading);
            if i == 0 {
                self.text_at(LEFT, y, Font::Regular, size, MUTED, &encode(label));
            }
            self.text_at(LEFT + indent, y, font, size, TEXT, line);
        }
    }

    fn gap(&mut self, h: f32) {
        self.y -= h;
    }

    /// The finished pages, each with its footer.
    fn finish(mut self, footer: &str) -> Vec<String> {
        if !self.page.is_empty() || self.pages.is_empty() {
            self.new_page();
        }
        let n = self.pages.len();
        let footer = encode(footer);
        let mut pages = std::mem::take(&mut self.pages);
        for (i, page) in pages.iter_mut().enumerate() {
            let mut d = Doc {
                pages: Vec::new(),
                page: std::mem::take(page),
                y: 0.0,
            };
            d.rule(BOTTOM - 22.0);
            d.text_at(LEFT, BOTTOM - 34.0, Font::Regular, 8.0, MUTED, &footer);
            let num = encode(&format!("Page {} of {}", i + 1, n));
            let w = width(Font::Regular, 8.0, &num);
            d.text_at(RIGHT - w, BOTTOM - 34.0, Font::Regular, 8.0, MUTED, &num);
            *page = d.page;
        }
        pages
    }
}

/// The line in the model's text a finding belongs to, as on the page.
fn line_of(r: &Report, f: &Finding) -> Option<usize> {
    f.line.or_else(|| {
        r.antimony
            .as_ref()
            .and_then(|a| a.line_of(f.element.kind(), &f.element.id()))
    })
}

fn plural(n: usize, one: &str, many: &str) -> String {
    format!("{n} {}", if n == 1 { one } else { many })
}

/// The report as a PDF file.
pub fn report(r: &Report, meta: &Meta) -> Vec<u8> {
    let name = r
        .model_name
        .clone()
        .or_else(|| r.model_id.clone())
        .unwrap_or_else(|| meta.file_name.to_string());
    let mut d = Doc::new();

    // ---- Title and what was checked
    let y = d.next_line(22.0);
    d.text_at(
        LEFT,
        y,
        Font::Bold,
        20.0,
        TEXT,
        &encode("Model check report"),
    );
    d.gap(8.0);
    // A model with neither name nor id is named by its file, shown below.
    let mut model = if r.model_name.is_some() || r.model_id.is_some() {
        name.clone()
    } else {
        String::new()
    };
    if let Some(id) = &r.model_id
        && *id != name
    {
        model.push_str(&format!(" (id {id})"));
    }
    for (label, value) in [
        ("Model", model.as_str()),
        ("File", meta.file_name),
        ("Format", meta.source),
        ("Created", meta.date),
        ("Checker", meta.checker),
    ] {
        if !value.is_empty() {
            d.labelled(label, 60.0, Font::Regular, 10.0, 14.0, value);
        }
    }

    // ---- Summary
    d.gap(10.0);
    let y = d.next_line(16.0);
    d.text_at(LEFT, y, Font::Bold, 13.0, TEXT, &encode("Summary"));
    d.gap(2.0);
    let counts = [
        (
            Severity::Error,
            plural(r.count(Severity::Error), "error", "errors"),
            "the model fails: a simulation will stop or give wrong results",
        ),
        (
            Severity::Warning,
            plural(r.count(Severity::Warning), "warning", "warnings"),
            "probably wrong, or could not be shown to be right",
        ),
        (
            Severity::Info,
            plural(
                r.count(Severity::Info),
                "information item",
                "information items",
            ),
            "worth knowing; often intended",
        ),
    ];
    for (sev, count, meaning) in &counts {
        let y = d.next_line(14.0);
        let c = encode(count);
        d.text_at(LEFT, y, Font::Bold, 10.0, colour(*sev), &c);
        let x = LEFT + 130.0;
        d.text_at(x, y, Font::Regular, 9.5, MUTED, &encode(meaning));
    }
    if r.findings.is_empty() {
        d.gap(4.0);
        d.para(
            LEFT,
            Font::Regular,
            10.0,
            TEXT,
            14.0,
            "No findings. The model passed every check.",
        );
    }

    // ---- Findings, by severity
    for (sev, heading) in [
        (Severity::Error, "Errors"),
        (Severity::Warning, "Warnings"),
        (Severity::Info, "Information"),
    ] {
        let found: Vec<&Finding> = r.findings.iter().filter(|f| f.severity == sev).collect();
        if found.is_empty() {
            continue;
        }
        d.gap(14.0);
        d.need(60.0);
        let y = d.next_line(16.0);
        d.text_at(
            LEFT,
            y,
            Font::Bold,
            13.0,
            colour(sev),
            &encode(&format!("{heading} ({})", found.len())),
        );
        d.gap(4.0);
        for f in found {
            finding(&mut d, r, f);
        }
    }

    // ---- The model, with line numbers
    if let Some(ant) = r.antimony.as_ref().filter(|a| !a.text.trim().is_empty()) {
        // The worst finding on each line, to colour its number.
        let mut worst: HashMap<usize, Severity> = HashMap::new();
        for f in &r.findings {
            if let Some(l) = line_of(r, f) {
                let e = worst.entry(l).or_insert(f.severity);
                if f.severity < *e {
                    *e = f.severity;
                }
            }
        }
        d.fresh_page();
        let y = d.next_line(16.0);
        d.text_at(
            LEFT,
            y,
            Font::Bold,
            13.0,
            TEXT,
            &encode("The model in Antimony"),
        );
        d.gap(2.0);
        d.para(
            LEFT,
            Font::Regular,
            9.0,
            MUTED,
            12.0,
            "Line numbers in colour have findings; the colour is the most serious one.",
        );
        d.gap(6.0);
        let size = 8.0;
        let indent = 34.0;
        for (i, text) in ant.text.lines().enumerate() {
            let n = i + 1;
            let parts = wrap(Font::Mono, size, RIGHT - LEFT - indent, &encode(text));
            for (j, part) in parts.iter().enumerate() {
                let y = d.next_line(10.5);
                if j == 0 {
                    let num = encode(&n.to_string());
                    let w = width(Font::Mono, size, &num);
                    let (font, c) = match worst.get(&n) {
                        Some(s) => (Font::Bold, colour(*s)),
                        None => (Font::Mono, MUTED),
                    };
                    d.text_at(LEFT + indent - 8.0 - w, y, font, size, c, &num);
                }
                d.text_at(LEFT + indent, y, Font::Mono, size, TEXT, part);
            }
        }
    }

    let footer = format!("Model check report: {name}");
    let pages = d.finish(&footer);
    let title = format!("Model check report: {name}");
    assemble(&pages, &title, meta.checker)
}

/// One finding: severity, code, element and line; then its title, the
/// reaction, the message, the expression and example values.
fn finding(d: &mut Doc, r: &Report, f: &Finding) {
    // Keep at least the head and a line of text together.
    d.need(46.0);
    let sev = f.severity;
    let y = d.next_line(13.0);
    let label = encode(&sev.as_str().to_uppercase());
    d.text_at(LEFT, y, Font::Bold, 8.5, colour(sev), &label);
    let mut x = LEFT + width(Font::Bold, 8.5, &label) + 8.0;
    let code = encode(f.code);
    d.text_at(x, y, Font::Mono, 8.5, MUTED, &code);
    x += width(Font::Mono, 8.5, &code) + 8.0;
    let mut place = f.element.to_string();
    if let Some(l) = line_of(r, f) {
        place.push_str(&format!(", line {l}"));
    }
    d.text_at(x, y, Font::Regular, 8.5, MUTED, &encode(&place));

    let title = f
        .title
        .as_deref()
        .unwrap_or_else(|| codes::lookup(f.code).map(|c| c.title).unwrap_or(""));
    if !title.is_empty() {
        d.para(LEFT, Font::Bold, 10.5, TEXT, 14.0, title);
    }

    // The message's last line is the expression, when it has one; it is
    // shown apart, and not at all when it repeats the reaction's rate law.
    let mut lines: Vec<&str> = f.message.split('\n').collect();
    let mut expr = None;
    if lines.len() > 1 {
        let last = lines[lines.len() - 1];
        if let Some(e) = last
            .strip_prefix("Rate law: ")
            .or_else(|| last.strip_prefix("Expression: "))
        {
            expr = Some(e);
            lines.pop();
        }
    }
    let law = f
        .reaction
        .as_deref()
        .and_then(|rx| rx.split_once("; ").map(|(_, l)| l));
    if expr.is_some() && expr == law {
        expr = None;
    }

    if let Some(rx) = &f.reaction {
        d.labelled("Reaction", 62.0, Font::Mono, 8.5, 11.5, rx);
    }
    d.gap(1.0);
    d.para(LEFT, Font::Regular, 9.5, TEXT, 12.5, &lines.join("\n"));
    if let Some(e) = expr {
        d.labelled("Expression", 62.0, Font::Mono, 8.5, 11.5, e);
    }
    if let Some(w) = &f.witness {
        d.labelled("Example", 62.0, Font::Mono, 8.5, 11.5, w);
    }
    d.gap(6.0);
    let y = d.y;
    d.rule(y);
    d.gap(4.0);
}

/// Puts the pages' content streams into a PDF file.
fn assemble(pages: &[String], title: &str, producer: &str) -> Vec<u8> {
    // Objects: 1 catalog, 2 page tree, 3-5 fonts, 6 info, then a page and
    // its content for each page.
    let n = pages.len();
    let page_id = |i: usize| 7 + 2 * i;
    let mut objects: Vec<Vec<u8>> = Vec::new();
    objects.push(b"<< /Type /Catalog /Pages 2 0 R >>".to_vec());
    let kids: Vec<String> = (0..n).map(|i| format!("{} 0 R", page_id(i))).collect();
    objects.push(format!("<< /Type /Pages /Kids [{}] /Count {n} >>", kids.join(" ")).into_bytes());
    for base in ["Helvetica", "Helvetica-Bold", "Courier"] {
        objects.push(
            format!(
                "<< /Type /Font /Subtype /Type1 /BaseFont /{base} /Encoding /WinAnsiEncoding >>"
            )
            .into_bytes(),
        );
    }
    objects.push(
        format!(
            "<< /Title {} /Producer {} /Creator {} >>",
            hex(&encode(title)),
            hex(&encode(producer)),
            hex(&encode(producer))
        )
        .into_bytes(),
    );
    for (i, content) in pages.iter().enumerate() {
        objects.push(
            format!(
                "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {PAGE_W} {PAGE_H}] \
                 /Resources << /Font << /F1 3 0 R /F2 4 0 R /F3 5 0 R >> >> \
                 /Contents {} 0 R >>",
                page_id(i) + 1
            )
            .into_bytes(),
        );
        let mut stream = format!("<< /Length {} >>\nstream\n", content.len()).into_bytes();
        stream.extend_from_slice(content.as_bytes());
        stream.extend_from_slice(b"\nendstream");
        objects.push(stream);
    }

    let mut out: Vec<u8> = b"%PDF-1.4\n%\xE2\xE3\xCF\xD3\n".to_vec();
    let mut offsets = Vec::with_capacity(objects.len());
    for (i, obj) in objects.iter().enumerate() {
        offsets.push(out.len());
        out.extend_from_slice(format!("{} 0 obj\n", i + 1).as_bytes());
        out.extend_from_slice(obj);
        out.extend_from_slice(b"\nendobj\n");
    }
    let xref = out.len();
    out.extend_from_slice(format!("xref\n0 {}\n", objects.len() + 1).as_bytes());
    out.extend_from_slice(b"0000000000 65535 f \n");
    for o in offsets {
        out.extend_from_slice(format!("{o:010} 00000 n \n").as_bytes());
    }
    out.extend_from_slice(
        format!(
            "trailer\n<< /Size {} /Root 1 0 R /Info 6 0 R >>\nstartxref\n{xref}\n%%EOF\n",
            objects.len() + 1
        )
        .as_bytes(),
    );
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Text as it appears in a page's content: hex, without the brackets.
    fn shown(s: &str) -> String {
        let h = hex(s.as_bytes());
        h[1..h.len() - 1].to_string()
    }

    fn meta() -> Meta<'static> {
        Meta {
            file_name: "model.xml",
            source: "SBML Level 3 Version 2",
            date: "17 Nov 2026",
            checker: "Model Checker test",
        }
    }

    #[test]
    fn wraps_at_spaces_and_breaks_long_words() {
        let lines = wrap(Font::Mono, 10.0, 60.0, b"aaaa bbbb cccccccccccccc");
        // 10 characters of Courier at 10 pt fit in 60 pt.
        assert_eq!(
            lines,
            vec![
                b"aaaa bbbb".to_vec(),
                b"cccccccccc".to_vec(),
                b"cccc".to_vec()
            ]
        );
        assert!(lines.iter().all(|l| width(Font::Mono, 10.0, l) <= 60.0));
    }

    #[test]
    fn encodes_as_win_ansi() {
        assert_eq!(
            encode("a\u{2192}b \u{b5}M \u{2013}"),
            b"a->b \xb5M \x96".to_vec()
        );
    }

    #[test]
    fn writes_a_well_formed_pdf() {
        let xml = r#"<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
  <model id="m">
    <listOfCompartments><compartment id="c" size="1" constant="true"/></listOfCompartments>
    <listOfSpecies>
      <species id="S" compartment="c" initialConcentration="1" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
      <species id="P" compartment="c" initialConcentration="0" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
    </listOfSpecies>
    <listOfParameters><parameter id="k" value="1" constant="true"/></listOfParameters>
    <listOfReactions>
      <reaction id="J1" reversible="false">
        <listOfReactants><speciesReference species="S" stoichiometry="1" constant="true"/></listOfReactants>
        <listOfProducts><speciesReference species="P" stoichiometry="1" constant="true"/></listOfProducts>
        <kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><apply><divide/><ci>k</ci><ci>P</ci></apply></math></kineticLaw>
      </reaction>
    </listOfReactions>
  </model>
</sbml>"#;
        let r = crate::check(xml);
        assert!(r.count(Severity::Error) > 0);
        let pdf = report(&r, &meta());
        assert!(pdf.starts_with(b"%PDF-1.4"));
        assert!(pdf.ends_with(b"%%EOF\n"));
        // Every cross-reference offset points at its object.
        let text = String::from_utf8_lossy(&pdf);
        let start: usize = text
            .rsplit("startxref\n")
            .next()
            .unwrap()
            .lines()
            .next()
            .unwrap()
            .parse()
            .unwrap();
        let tail = String::from_utf8_lossy(&pdf[start..]).into_owned();
        let entries: Vec<&str> = tail.lines().skip(3).collect();
        let mut n = 1;
        for e in entries {
            let Some(off) = e.strip_suffix(" 00000 n ") else {
                break;
            };
            let off: usize = off.parse().unwrap();
            assert!(pdf[off..].starts_with(format!("{n} 0 obj").as_bytes()));
            n += 1;
        }
        assert!(n > 7);
        // The finding's code is in the text (hex encoded).
        assert!(text.contains(&shown("RL012")));
    }

    #[test]
    fn many_findings_run_onto_more_pages() {
        let mut r = crate::check("<not sbml");
        let f = r.findings[0].clone();
        r.findings = vec![f; 120];
        let pdf = report(&r, &meta());
        let text = String::from_utf8_lossy(&pdf);
        assert!(text.contains(&shown("Page 2 of")));
    }
}
