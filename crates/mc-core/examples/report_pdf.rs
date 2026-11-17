//! Writes the PDF report the web page saves, for an SBML or Antimony file.
//!
//! cargo run -p mc-core --features antimony --example report_pdf -- model.xml > report.pdf

use mc_core::{from_antimony, pdf};

fn main() {
    let Some(path) = std::env::args().nth(1) else {
        eprintln!("usage: report_pdf FILE (SBML or Antimony)");
        std::process::exit(2);
    };
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| {
        eprintln!("{path}: {e}");
        std::process::exit(1);
    });
    // As on the page: SBML is XML, so starts with "<"; Antimony never does.
    let antimony = !text
        .trim_start_matches('\u{feff}')
        .trim_start()
        .starts_with('<');
    let report = if antimony {
        from_antimony::check_antimony(&text)
    } else {
        mc_core::check(&text)
    };
    let source = match (antimony, report.level, report.version) {
        (true, _, _) => "Antimony".to_string(),
        (false, Some(l), Some(v)) => format!("SBML Level {l} Version {v}"),
        _ => "SBML".to_string(),
    };
    let file_name = std::path::Path::new(&path)
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let checker = format!("Model Checker {}", env!("CARGO_PKG_VERSION"));
    let bytes = pdf::report(
        &report,
        &pdf::Meta {
            file_name: &file_name,
            source: &source,
            date: &now_utc(),
            checker: &checker,
        },
    );
    use std::io::Write;
    std::io::stdout().write_all(&bytes).unwrap();
}

/// The current date and time in UTC, as "17 November 2026, 21:15 UTC".
fn now_utc() -> String {
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    let (days, rest) = (secs.div_euclid(86400), secs.rem_euclid(86400));
    // Days since 1970-01-01 to a civil date (Howard Hinnant's algorithm).
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = if mp < 10 { mp + 3 } else { mp - 9 };
    let year = yoe + era * 400 + i64::from(month <= 2);
    const MONTHS: [&str; 12] = [
        "January",
        "February",
        "March",
        "April",
        "May",
        "June",
        "July",
        "August",
        "September",
        "October",
        "November",
        "December",
    ];
    format!(
        "{day} {} {year}, {:02}:{:02} UTC",
        MONTHS[(month - 1) as usize],
        rest / 3600,
        rest % 3600 / 60
    )
}
