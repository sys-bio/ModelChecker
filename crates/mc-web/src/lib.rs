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
