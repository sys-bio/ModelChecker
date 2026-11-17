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
