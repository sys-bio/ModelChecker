//! The checks, one module per family.

mod annotations;
mod features;
mod mass_action;
mod network;
mod poly;
mod reactions;
mod rules;
mod sbo;
pub mod scope;
mod singular;
mod symbols;
mod thermo;
mod units;
mod values;

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
            reaction: None,
            title: None,
            line: None,
        });
    }

    /// As [`Ctx::add`], with a title for this particular problem.
    pub fn add_titled(
        &mut self,
        code: &'static str,
        severity: Severity,
        element: Element,
        title: String,
        message: impl Into<String>,
        witness: Option<String>,
    ) {
        self.add(code, severity, element, message, witness);
        if let Some(f) = self.findings.last_mut() {
            f.title = Some(title);
        }
    }
}

/// Runs every check on `model`. `namespaces` are those declared on the
/// document's `<sbml>` element, which annotations may use undeclared.
pub fn run(model: &Model, level: u32, namespaces: &[(Option<String>, String)]) -> Vec<Finding> {
    let view = View::new(model, level);
    let mut ctx = Ctx {
        view: &view,
        findings: Vec::new(),
    };
    annotations::check(&mut ctx, namespaces);
    features::check(&mut ctx, namespaces);
    symbols::check(&mut ctx);
    singular::check(&mut ctx);
    reactions::check(&mut ctx);
    mass_action::check(&mut ctx);
    thermo::check(&mut ctx);
    network::check(&mut ctx);
    rules::check(&mut ctx);
    values::check(&mut ctx);
    units::check(&mut ctx);
    sbo::check(&mut ctx);
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
