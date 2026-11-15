//! The checker for the browser. The page calls these through
//! `window.wasmBindings`, which Trunk sets up.

use wasm_bindgen::prelude::*;

/// Checks an SBML document; returns the report as JSON.
#[wasm_bindgen]
pub fn check_sbml(xml: &str) -> String {
    mc_core::json::report(&mc_core::check(xml))
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
