//! The checks, one module per family.

mod reactions;
pub mod scope;
mod singular;
mod symbols;
mod units;

use std::collections::{HashMap, HashSet};

use sbml_io::{Math, Model, Reaction};

use crate::finding::{Element, Finding, Severity};
use crate::view::View;

/// Where checks put what they find.
pub struct Ctx<'v, 'm> {
    pub view: &'v View<'m>,
    pub findings: Vec<Finding>,
}

impl Ctx<'_, '_> {
    pub fn add(
        &mut self,
        code: &'static str,
        severity: Severity,
        element: Element,
        message: impl Into<String>,
        witness: Option<String>,
    ) {
        self.findings.push(Finding {
            code,
            severity,
            element,
            message: message.into(),
            witness,
        });
    }
}

/// Runs every check on `model`.
pub fn run(model: &Model, level: u32) -> Vec<Finding> {
    let view = View::new(model, level);
    let mut ctx = Ctx {
        view: &view,
        findings: Vec::new(),
    };
    symbols::check(&mut ctx);
    singular::check(&mut ctx);
    reactions::check(&mut ctx);
    units::check(&mut ctx);
    ctx.findings
}

/// A reaction's local parameters: id to value (NaN if none).
pub fn locals(r: &Reaction) -> HashMap<String, f64> {
    r.kinetic_law
        .iter()
        .flat_map(|k| &k.parameters)
        .map(|p| (p.id.clone(), p.value.unwrap_or(f64::NAN)))
        .collect()
}

pub fn shadow(locals: &HashMap<String, f64>) -> HashSet<&str> {
    locals.keys().map(String::as_str).collect()
}

/// A reaction's rate law with functions and assignment rules inlined.
pub fn expanded_rate(view: &View, r: &Reaction) -> Option<Math> {
    let math = r.kinetic_law.as_ref()?.math.as_ref()?;
    let locals = locals(r);
    Some(view.expand(math, &shadow(&locals)))
}
