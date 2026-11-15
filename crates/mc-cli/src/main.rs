//! modelcheck: run the checks on SBML files or folders of them.
//!
//! ```text
//! modelcheck model.xml                 findings for one model
//! modelcheck --summary corpus/         counts per code over a folder
//! modelcheck --code RL003 corpus/      only findings with this code
//! modelcheck --json model.xml          the report as JSON
//! ```

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::Instant;

use mc_core::{Report, Severity, codes};

struct Options {
    summary: bool,
    json: bool,
    code: Option<String>,
    quiet_info: bool,
    paths: Vec<PathBuf>,
}

fn usage() -> ExitCode {
    eprintln!(
        "usage: modelcheck [--summary] [--json] [--no-info] [--code CODE] FILE_OR_DIR...\n\
         \n  --summary   counts per code instead of each finding\
         \n  --json      the report as JSON (single file)\
         \n  --no-info   leave out info findings\
         \n  --code C    only findings with code C"
    );
    ExitCode::from(2)
}

fn main() -> ExitCode {
    let mut o = Options {
        summary: false,
        json: false,
        code: None,
        quiet_info: false,
        paths: Vec::new(),
    };
    let mut args = std::env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--summary" => o.summary = true,
            "--json" => o.json = true,
            "--no-info" => o.quiet_info = true,
            "--code" => match args.next() {
                Some(c) => o.code = Some(c),
                None => return usage(),
            },
            "-h" | "--help" => return usage(),
            _ => o.paths.push(PathBuf::from(a)),
        }
    }
    if o.paths.is_empty() {
        return usage();
    }
    let mut files = Vec::new();
    for p in &o.paths {
        collect(p, &mut files);
    }
    files.sort();

    // code -> (findings, files with at least one)
    let mut totals: BTreeMap<&'static str, (usize, usize)> = BTreeMap::new();
    let mut errors = 0usize;
    let start = Instant::now();
    for f in &files {
        let text = match std::fs::read_to_string(f) {
            Ok(t) => t,
            Err(e) => {
                eprintln!("{}: cannot read: {e}", f.display());
                continue;
            }
        };
        let report = mc_core::check(&text);
        if o.json {
            println!("{}", mc_core::json::report(&report));
            continue;
        }
        let shown: Vec<_> = report
            .findings
            .iter()
            .filter(|x| !(o.quiet_info && x.severity == Severity::Info))
            .filter(|x| o.code.as_deref().is_none_or(|c| c == x.code))
            .collect();
        errors += shown
            .iter()
            .filter(|x| x.severity == Severity::Error)
            .count();
        let mut seen_codes = Vec::new();
        for x in &shown {
            let t = totals.entry(x.code).or_default();
            t.0 += 1;
            if !seen_codes.contains(&x.code) {
                seen_codes.push(x.code);
                t.1 += 1;
            }
        }
        if !o.summary && !shown.is_empty() {
            print_report(f, &report, &shown);
        }
    }
    if o.summary || files.len() > 1 {
        println!(
            "\n{} files in {:.1} s",
            files.len(),
            start.elapsed().as_secs_f64()
        );
        println!("{:<7} {:>8} {:>7}  title", "code", "findings", "files");
        for (code, (n, nf)) in &totals {
            let title = codes::lookup(code).map(|c| c.title).unwrap_or("");
            println!("{code:<7} {n:>8} {nf:>7}  {title}");
        }
    }
    if errors > 0 {
        ExitCode::from(1)
    } else {
        ExitCode::SUCCESS
    }
}

fn print_report(path: &Path, r: &Report, shown: &[&mc_core::Finding]) {
    let name = r
        .model_id
        .as_deref()
        .map(|id| format!(" ({id})"))
        .unwrap_or_default();
    println!("== {}{name}", path.display());
    for f in shown {
        let mut lines = f.message.lines();
        println!(
            "  {:<7} {} [{}]: {}",
            f.severity.as_str(),
            f.code,
            f.element,
            lines.next().unwrap_or("")
        );
        if let Some(r) = &f.reaction {
            println!("          reaction: {r}");
        }
        // The rate law as written is already in the reaction line; an
        // expanded one (functions, rules inlined) is still worth showing.
        let law = f
            .reaction
            .as_deref()
            .and_then(|r| r.split_once("; "))
            .map(|(_, l)| l);
        for l in lines {
            let text = l
                .strip_prefix("Rate law: ")
                .or_else(|| l.strip_prefix("Expression: "));
            if law.is_some() && text == law {
                continue;
            }
            println!("          {l}");
        }
        if let Some(w) = &f.witness {
            println!("          witness: {w}");
        }
    }
}

fn collect(p: &Path, out: &mut Vec<PathBuf>) {
    if p.is_dir() {
        let Ok(rd) = std::fs::read_dir(p) else { return };
        for e in rd.flatten() {
            collect(&e.path(), out);
        }
    } else if p
        .extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("xml") || e.eq_ignore_ascii_case("sbml"))
    {
        out.push(p.to_path_buf());
    }
}
