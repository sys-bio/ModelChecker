//! The checker for the browser. The page calls these through
//! `window.wasmBindings`, which Trunk sets up.

use wasm_bindgen::prelude::*;

/// Checks an SBML document; returns the report as JSON.
#[wasm_bindgen]
pub fn check_sbml(xml: &str) -> String {
    mc_core::json::report(&mc_core::check(xml))
}

/// Checks the original model (`original_xml`) with an edit of its Antimony
/// text applied; returns the report as JSON. The edit is merged into the
/// original, so nothing the text does not show is lost.
#[wasm_bindgen]
pub fn check_edit(original_xml: &str, edited: &str) -> String {
    mc_core::json::report(&mc_core::from_antimony::check_edit(original_xml, edited))
}

/// The original model with the edit applied, as SBML; an empty string if
/// the edited text has errors.
#[wasm_bindgen]
pub fn edited_sbml(original_xml: &str, edited: &str) -> String {
    mc_core::from_antimony::edited_sbml(original_xml, edited).unwrap_or_default()
}

/// Checks a model written in Antimony; returns the report as JSON. The
/// text is the master copy, so findings refer to its own lines.
#[wasm_bindgen]
pub fn check_antimony(text: &str) -> String {
    mc_core::json::report(&mc_core::from_antimony::check_antimony(text))
}

/// An Antimony model as SBML; an empty string if the text has errors.
#[wasm_bindgen]
pub fn antimony_sbml(text: &str) -> String {
    mc_core::from_antimony::antimony_sbml(text).unwrap_or_default()
}

/// The report of a check as a PDF file. `kind` says what was checked, as
/// for the functions above: "sbml" (`xml`), "edit" (`xml` with the edited
/// `text`) or "antimony" (`text`). `date` is when, as the reader should
/// see it.
#[wasm_bindgen]
pub fn report_pdf(kind: &str, xml: &str, text: &str, file_name: &str, date: &str) -> Vec<u8> {
    use mc_core::from_antimony;
    let report = match kind {
        "antimony" => from_antimony::check_antimony(text),
        "edit" => from_antimony::check_edit(xml, text),
        _ => mc_core::check(xml),
    };
    let sbml = match (report.level, report.version) {
        (Some(l), Some(v)) => format!("SBML Level {l} Version {v}"),
        _ => "SBML".to_string(),
    };
    let source = match kind {
        "antimony" => "Antimony".to_string(),
        "edit" => format!("{sbml}, edited in the Antimony panel"),
        _ => sbml,
    };
    let commit = env!("MC_BUILD_COMMIT");
    let checker = format!(
        "Model Checker {}{}, https://sys-bio.github.io/ModelChecker/",
        env!("CARGO_PKG_VERSION"),
        if commit.len() >= 7 {
            format!(" (build {})", &commit[..7])
        } else {
            String::new()
        }
    );
    mc_core::pdf::report(
        &report,
        &mc_core::pdf::Meta {
            file_name,
            source: &source,
            date,
            checker: &checker,
        },
    )
}

/// Every finding code with its title and explanation, as JSON.
#[wasm_bindgen]
pub fn finding_codes() -> String {
    mc_core::json::codes()
}

/// The checker's version.
#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_string()
}

/// The commit the page was built from (full hash; empty if unknown).
#[wasm_bindgen]
pub fn build_commit() -> String {
    env!("MC_BUILD_COMMIT").to_string()
}

/// The date the page was built, as YYYY-MM-DD (UTC).
#[wasm_bindgen]
pub fn build_date() -> String {
    env!("MC_BUILD_DATE").to_string()
}

/// Whether the page was built by GitHub Actions (the published site)
/// rather than locally.
#[wasm_bindgen]
pub fn built_by_ci() -> bool {
    env!("MC_BUILD_CI") == "true"
}
