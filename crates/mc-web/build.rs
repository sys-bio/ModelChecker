//! Records which commit the page was built from, and when, for the footer.
//!
//! On GitHub Actions the commit is `GITHUB_SHA`; elsewhere it is asked of
//! git, and the build is marked local. The date is the build's (UTC).

use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn main() {
    let ci = std::env::var("GITHUB_ACTIONS").as_deref() == Ok("true");
    let commit = std::env::var("GITHUB_SHA")
        .ok()
        .filter(|s| !s.is_empty())
        .or_else(|| {
            Command::new("git")
                .args(["rev-parse", "HEAD"])
                .output()
                .ok()
                .filter(|o| o.status.success())
                .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
        })
        .unwrap_or_default();
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    let (y, m, d) = civil_from_days((secs / 86_400) as i64);
    println!("cargo:rustc-env=MC_BUILD_COMMIT={commit}");
    println!("cargo:rustc-env=MC_BUILD_DATE={y:04}-{m:02}-{d:02}");
    println!("cargo:rustc-env=MC_BUILD_CI={ci}");
    // Rebuild when the commit changes (a commit, a checkout) or in CI.
    println!("cargo:rerun-if-env-changed=GITHUB_SHA");
    println!("cargo:rerun-if-changed=../../.git/HEAD");
    println!("cargo:rerun-if-changed=../../.git/logs/HEAD");
}

/// Days since 1970-01-01 to (year, month, day): Howard Hinnant's algorithm.
fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let y = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    (if m <= 2 { y + 1 } else { y }, m, d)
}
