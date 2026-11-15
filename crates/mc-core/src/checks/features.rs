//! Features many simulators do not support, or handle differently.
//!
//! - FT001 `fast` reactions (most simulators treat them as ordinary ones;
//!   SBML Level 3 Version 2 dropped the attribute);
//! - FT002 delays (`delay(x, t)`): a delay differential equation;
//! - FT003 newer mathematics: `rateOf`, and `min`, `max`, `rem`,
//!   `quotient` and `implies` (Level 3 Version 2);
//! - FT004 SBML packages: `comp` (hierarchical models; this checker sees
//!   only the top model) and others that change what the model means;
//! - FT005 variable stoichiometry;
//! - FT006 constraints.

use std::collections::BTreeSet;

use sbml_io::{Math, Op, Symbol};

use super::Ctx;
use crate::codes::*;
use crate::finding::{Element, Severity};

pub fn check(ctx: &mut Ctx, namespaces: &[(Option<String>, String)]) {
    fast(ctx);
    maths(ctx);
    packages(ctx, namespaces);
    variable_stoichiometry(ctx);
    constraints(ctx);
}

fn list(items: &[String], n: usize) -> String {
    if items.len() <= n {
        items.join(", ")
    } else {
        format!("{}, ... ({} in all)", items[..n].join(", "), items.len())
    }
}

/// FT001.
fn fast(ctx: &mut Ctx) {
    let fast: Vec<String> = ctx
        .view
        .model
        .reactions
        .iter()
        .filter(|r| r.fast == Some(true))
        .map(|r| r.id.clone())
        .collect();
    if fast.is_empty() {
        return;
    }
    ctx.add(
        FT_FAST,
        Severity::Warning,
        Element::Model,
        format!(
            "{} marked fast ({}). Most simulators ignore the mark and treat {}, so results \
             differ between tools; SBML Level 3 Version 2 removed it. To be safe, write the \
             fast step as a rapid equilibrium explicitly (for example with an assignment rule \
             or large rate constants).",
            if fast.len() == 1 { "One reaction is".to_string() } else { format!("{} reactions are", fast.len()) },
            list(&fast, 10),
            if fast.len() == 1 { "it as an ordinary reaction" } else { "them as ordinary reactions" }
        ),
        None,
    );
}

/// Every expression in the model, with the element it belongs to.
fn expressions(m: &sbml_io::Model) -> Vec<(Element, &Math)> {
    let mut out: Vec<(Element, &Math)> = Vec::new();
    for f in &m.function_definitions {
        if let Some(x) = &f.math {
            out.push((Element::FunctionDefinition(f.id.clone()), x));
        }
    }
    for r in &m.reactions {
        if let Some(x) = r.kinetic_law.as_ref().and_then(|k| k.math.as_ref()) {
            out.push((Element::Reaction(r.id.clone()), x));
        }
    }
    for (i, rule) in m.rules.iter().enumerate() {
        if let Some(x) = rule.math() {
            let el = match rule.variable() {
                Some(v) => Element::Rule(v.to_string()),
                None => Element::AlgebraicRule(i),
            };
            out.push((el, x));
        }
    }
    for ia in &m.initial_assignments {
        if let Some(x) = &ia.math {
            out.push((Element::InitialAssignment(ia.symbol.clone()), x));
        }
    }
    for (i, e) in m.events.iter().enumerate() {
        let el = Element::Event(e.id.clone().unwrap_or_else(|| format!("#{}", i + 1)));
        for x in e
            .trigger
            .iter()
            .filter_map(|t| t.math.as_ref())
            .chain(e.delay.iter().filter_map(|d| d.math.as_ref()))
            .chain(e.priority.iter().filter_map(|p| p.math.as_ref()))
            .chain(e.assignments.iter().filter_map(|a| a.math.as_ref()))
        {
            out.push((el.clone(), x));
        }
    }
    out
}

/// The delays, `rateOf`s and Level 3 Version 2 operators in `m`.
fn visit(m: &Math, delay: &mut bool, rate_of: &mut bool, newer: &mut BTreeSet<&'static str>) {
    match m {
        Math::CallSymbol(s, args) => {
            if s.url == Symbol::DELAY {
                *delay = true;
            } else if s.url == Symbol::RATE_OF {
                *rate_of = true;
            }
            args.iter().for_each(|a| visit(a, delay, rate_of, newer));
        }
        Math::Apply(op, args) => {
            match op {
                Op::Min => {
                    newer.insert("min");
                }
                Op::Max => {
                    newer.insert("max");
                }
                Op::Rem => {
                    newer.insert("rem");
                }
                Op::Quotient => {
                    newer.insert("quotient");
                }
                Op::Implies => {
                    newer.insert("implies");
                }
                _ => {}
            }
            args.iter().for_each(|a| visit(a, delay, rate_of, newer));
        }
        Math::Call(_, args) => args.iter().for_each(|a| visit(a, delay, rate_of, newer)),
        Math::Root { degree, arg } | Math::Log { base: degree, arg } => {
            if let Some(d) = degree {
                visit(d, delay, rate_of, newer);
            }
            visit(arg, delay, rate_of, newer);
        }
        Math::Piecewise { pieces, otherwise } => {
            for (v, c) in pieces {
                visit(v, delay, rate_of, newer);
                visit(c, delay, rate_of, newer);
            }
            if let Some(o) = otherwise {
                visit(o, delay, rate_of, newer);
            }
        }
        Math::Lambda { body, .. } => visit(body, delay, rate_of, newer),
        Math::Semantics { math, .. } => visit(math, delay, rate_of, newer),
        _ => {}
    }
}

/// FT002 and FT003.
fn maths(ctx: &mut Ctx) {
    let mut delayed: Vec<String> = Vec::new();
    let mut rate_of: Vec<String> = Vec::new();
    let mut newer: BTreeSet<&'static str> = BTreeSet::new();
    let mut newer_where: Vec<String> = Vec::new();
    let model = ctx.view.model;
    for (el, m) in expressions(model) {
        let (mut d, mut r) = (false, false);
        let mut n = BTreeSet::new();
        visit(m, &mut d, &mut r, &mut n);
        let name = el.to_string();
        if d && !delayed.contains(&name) {
            delayed.push(name.clone());
        }
        if r && !rate_of.contains(&name) {
            rate_of.push(name.clone());
        }
        if !n.is_empty() && !newer_where.contains(&name) {
            newer_where.push(name);
        }
        newer.extend(n);
    }
    if !delayed.is_empty() {
        ctx.add(
            FT_DELAY,
            Severity::Info,
            Element::Model,
            format!(
                "The model uses delays (`delay(x, t)`, a value from some time ago) in {}. That \
                 makes it a delay differential equation, which many simulators cannot solve.",
                list(&delayed, 8)
            ),
            None,
        );
    }
    let mut parts: Vec<String> = Vec::new();
    if !rate_of.is_empty() {
        parts.push(format!("`rateOf` (in {})", list(&rate_of, 6)));
    }
    if !newer.is_empty() {
        let ops: Vec<String> = newer.iter().map(|s| format!("`{s}`")).collect();
        parts.push(format!("{} (in {})", ops.join(", "), list(&newer_where, 6)));
    }
    if !parts.is_empty() {
        ctx.add(
            FT_NEWER_MATH,
            Severity::Info,
            Element::Model,
            format!(
                "The model uses mathematics added in SBML Level 3 Version 2: {}. Tools that \
                 read only older versions of SBML cannot read these.",
                parts.join("; ")
            ),
            None,
        );
    }
}

/// FT004.
fn packages(ctx: &mut Ctx, namespaces: &[(Option<String>, String)]) {
    // Packages that change what the model means. Layout, render, groups
    // and fbc do not change a simulation.
    const PACKAGES: &[(&str, &str, &str)] = &[
        ("/comp/", "comp", "hierarchical model composition (submodels)"),
        ("/distrib/", "distrib", "random distributions"),
        ("/arrays/", "arrays", "arrays of elements"),
        ("/multi/", "multi", "multistate, multicomponent species"),
        ("/spatial/", "spatial", "spatial geometry"),
        ("/qual/", "qual", "qualitative (logical) models"),
    ];
    let mut seen: BTreeSet<&str> = BTreeSet::new();
    for (_, uri) in namespaces {
        for (marker, name, what) in PACKAGES {
            if uri.contains(marker) && seen.insert(name) {
                let (severity, message) = if *name == "comp" {
                    (
                        Severity::Warning,
                        "The model uses the SBML comp package (hierarchical model composition): \
                         parts of it may live in submodels. This checker looks only at the top \
                         model, so findings may be missing or wrong, and many simulators need \
                         the model flattened first."
                            .to_string(),
                    )
                } else {
                    (
                        Severity::Info,
                        format!(
                            "The model uses the SBML {name} package ({what}). This checker does \
                             not analyse it, and many simulators do not support it."
                        ),
                    )
                };
                ctx.add(FT_PACKAGE, severity, Element::Model, message, None);
            }
        }
    }
}

/// FT005.
fn variable_stoichiometry(ctx: &mut Ctx) {
    let view = ctx.view;
    let set: BTreeSet<&str> = view
        .assignment_rules
        .keys()
        .copied()
        .chain(view.rate_rules.iter().copied())
        .chain(view.event_targets.iter().copied())
        .chain(view.initial_assignments.keys().copied())
        .collect();
    let mut variable: Vec<String> = Vec::new();
    for r in &view.model.reactions {
        let changes = r.reactants.iter().chain(&r.products).any(|sr| {
            sr.stoichiometry_math.is_some() || sr.id.as_deref().is_some_and(|id| set.contains(id))
        });
        if changes {
            variable.push(r.id.clone());
        }
    }
    if variable.is_empty() {
        return;
    }
    ctx.add(
        FT_VARIABLE_STOICHIOMETRY,
        Severity::Info,
        Element::Model,
        format!(
            "Stoichiometries are given by mathematics or set by rules or events in {}. Some \
             simulators do not support this, and several of this checker's checks assume fixed \
             stoichiometries.",
            list(&variable, 10)
        ),
        None,
    );
}

/// FT006.
fn constraints(ctx: &mut Ctx) {
    let n = ctx.view.model.constraints.len();
    if n == 0 {
        return;
    }
    ctx.add(
        FT_CONSTRAINT,
        Severity::Info,
        Element::Model,
        format!(
            "The model has {} (conditions that should hold during a simulation). Few \
             simulators check them, so a run may break one without warning.",
            if n == 1 { "a constraint".to_string() } else { format!("{n} constraints") }
        ),
        None,
    );
}
