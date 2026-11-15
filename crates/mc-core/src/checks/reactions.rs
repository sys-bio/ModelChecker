//! Checks on each reaction against its rate law:
//!
//! - RL008 no kinetic law, RL011 no reactants or products;
//! - RL005 zero order, RL006 a reactant missing from the rate law,
//!   RL007 a constant-rate source;
//! - RL003/RL004 positivity: with a species at zero, the reaction must not
//!   consume it. For net stoichiometry `n` and rate `v`, `n * v >= 0` must
//!   hold at that species = 0, for all non-negative values of the rest.
//!   The `reversible` attribute is ignored: the rate law alone decides.
//! - RL009 a species in the rate law that is not declared, RL010 a
//!   modifier the rate law does not use.

use std::collections::{BTreeMap, HashSet};

use sbml_io::{Math, Reaction};

use super::scope::{BUDGET, Scope, free_symbols};
use super::{Ctx, expanded_rate, locals};
use crate::codes::*;
use crate::finding::{Element, Severity};
use crate::fmt::infix;
use crate::sample;
use crate::view::{Kind, idents};

pub fn check(ctx: &mut Ctx) {
    let view = ctx.view;
    for r in &view.model.reactions {
        let el = Element::Reaction(r.id.clone());
        if r.reactants.is_empty() && r.products.is_empty() {
            ctx.add(
                RL_EMPTY_REACTION,
                Severity::Warning,
                el.clone(),
                "The reaction has no reactants and no products.",
                None,
            );
        }
        let Some(rate) = expanded_rate(view, r) else {
            ctx.add(
                RL_NO_KINETIC_LAW,
                Severity::Warning,
                el,
                "The reaction has no kinetic law.",
                None,
            );
            continue;
        };
        declared(ctx, r);
        order_and_positivity(ctx, r, &rate);
    }
}

/// Net stoichiometry of each species the reaction changes.
fn net_stoichiometry(r: &Reaction) -> BTreeMap<String, f64> {
    let mut net: BTreeMap<String, f64> = BTreeMap::new();
    for (refs, sign) in [(&r.reactants, -1.0), (&r.products, 1.0)] {
        for sr in refs {
            // Stoichiometry given as math (Level 2) is taken as 1: only its
            // sign matters here, and it is almost always positive.
            let st = sr.stoichiometry.unwrap_or(1.0);
            *net.entry(sr.species.clone()).or_default() += sign * st;
        }
    }
    net
}

fn is_held(ctx: &Ctx, species: &str) -> bool {
    ctx.view
        .model
        .species
        .iter()
        .find(|s| s.id == species)
        .is_some_and(|s| s.boundary_condition == Some(true) || s.constant == Some(true))
}

fn order_and_positivity(ctx: &mut Ctx, r: &Reaction, rate: &Math) {
    let view = ctx.view;
    let el = Element::Reaction(r.id.clone());
    let deps: HashSet<String> = free_symbols(rate).into_iter().collect();
    let species_deps: Vec<&String> = deps
        .iter()
        .filter(|d| view.kinds.get(d.as_str()) == Some(&Kind::Species))
        .collect();
    let net = net_stoichiometry(r);
    // Species this reaction really changes.
    let changed: Vec<(&String, f64)> = net
        .iter()
        .filter(|(s, n)| **n != 0.0 && !is_held(ctx, s))
        .map(|(s, n)| (s, *n))
        .collect();
    let consumed: Vec<&String> = changed.iter().filter(|(_, n)| *n < 0.0).map(|(s, _)| *s).collect();
    let missing: Vec<&String> = consumed.iter().copied().filter(|s| !deps.contains(*s)).collect();

    // Zero order, missing reactants, constant sources.
    if !consumed.is_empty() && missing.len() == consumed.len() {
        let names = list(&consumed);
        ctx.add(
            RL_ZERO_ORDER,
            Severity::Warning,
            el.clone(),
            format!(
                "The reaction consumes {names} but its rate does not depend on {}. \
                 This may be intended (a saturated enzyme), but it is often a typo, \
                 and {names} can go negative.\nRate law: {}",
                if consumed.len() == 1 { "it" } else { "any of them" },
                infix(rate)
            ),
            None,
        );
    } else {
        for s in &missing {
            ctx.add(
                RL_MISSING_REACTANT,
                Severity::Warning,
                el.clone(),
                format!(
                    "The reaction consumes {s} but its rate does not depend on {s}, \
                     so {s} can go negative.\nRate law: {}",
                    infix(rate)
                ),
                None,
            );
        }
    }
    if consumed.is_empty() && species_deps.is_empty() {
        ctx.add(
            RL_CONSTANT_SOURCE,
            Severity::Info,
            el.clone(),
            format!(
                "The rate depends on no species: a zero-order source.\nRate law: {}",
                infix(rate)
            ),
            None,
        );
    }

    // Positivity, for every changed species except reactants already
    // reported as missing from the rate law.
    let mut base = Scope::new(view);
    base.locals = locals(r);
    for (s, n) in &changed {
        if missing.contains(s) {
            continue;
        }
        // Not a concentration (negative initial value): no sign to keep.
        if view.range(s).is_some_and(|r| r.interval.lo < 0.0) {
            continue;
        }
        positivity(ctx, &base, r, rate, s, *n);
    }
}

fn positivity(ctx: &mut Ctx, base: &Scope, r: &Reaction, rate: &Math, s: &str, n: f64) {
    let mut scope = Scope {
        view: base.view,
        locals: base.locals.clone(),
        pins: base.pins.clone(),
    };
    scope.pins.insert(s.to_string(), 0.0);
    let role = if n < 0.0 { "reactant" } else { "product" };
    let needed = if n < 0.0 { "<= 0" } else { ">= 0" };

    // Proved by intervals?
    let (iv, _) = scope.interval(rate);
    let ok = if n < 0.0 { iv.hi <= 0.0 } else { iv.lo >= 0.0 };
    if ok {
        return;
    }
    // Look for a counterexample: first with every other variable between
    // zero and its initial value, where a model whose totals are written as
    // parameters (`ET - ES`) stays consistent; then anywhere.
    let vars = scope.vars(rate);
    let extras = scope.extras(rate);
    let typical = typical_box(&vars);
    let dir = if n < 0.0 { "consumes" } else { "consumes (runs backwards on)" };
    if let Some(w) = counterexample(&scope, rate, &typical, &extras, s, n) {
        ctx.add(
            RL_NEGATIVITY,
            Severity::Warning,
            Element::Reaction(r.id.clone()),
            format!(
                "{s} is a net {role}, so the rate must be {needed} when {s} = 0, but it \
                 is not: the reaction {dir} {s} when there is none, and {s} can go \
                 negative.\nRate law: {}",
                infix(rate)
            ),
            Some(w),
        );
        return;
    }
    let wide = counterexample(&scope, rate, &vars, &extras, s, n);
    // Where the rate law may itself be singular its sign means little, and
    // the singularity is reported on its own.
    if wide.is_none() {
        let (_, hazards) = scope.interval(rate);
        if hazards
            .iter()
            .any(|h| h.status != crate::value::Status::Safe && !h.guarded)
        {
            return;
        }
    }
    let message = if wide.is_some() {
        format!(
            "{s} is a net {role}, so the rate must be {needed} when {s} = 0. That fails \
             only when other variables exceed their initial values, which may not be \
             reachable (for example a part larger than its total).\nRate law: {}",
            infix(rate)
        )
    } else {
        format!(
            "{s} is a net {role}, so the rate must be {needed} when {s} = 0. This could not \
             be proved, and no counterexample was found.\nRate law: {}",
            infix(rate)
        )
    };
    ctx.add(
        RL_POSSIBLE_NEGATIVITY,
        Severity::Info,
        Element::Reaction(r.id.clone()),
        message,
        wide,
    );
}

/// `vars` limited to between their lowest value and their initial value.
/// Time, which starts at zero, is allowed a little way forward.
fn typical_box(vars: &[sample::Var]) -> Vec<sample::Var> {
    vars.iter()
        .map(|v| {
            let mut v = v.clone();
            let lo = if v.interval.lo.is_finite() {
                v.interval.lo
            } else {
                v.nominal.min(0.0)
            };
            let hi = if v.id == crate::eval::TIME {
                100.0
            } else {
                v.nominal.max(lo)
            };
            v.interval = crate::value::Interval::new(lo, hi);
            v
        })
        .collect()
}

/// A reachable point, with `s` at zero, where `n * rate < 0`, as text.
fn counterexample(
    scope: &Scope,
    rate: &Math,
    vars: &[sample::Var],
    extras: &[f64],
    s: &str,
    n: f64,
) -> Option<String> {
    for p in sample::points(vars, extras, BUDGET) {
        if !scope.reachable(vars, &p) {
            continue;
        }
        let (v, _) = scope.at(rate, vars, &p);
        // A small absolute floor keeps round-off from counting.
        if v.is_finite() && n * v < 0.0 && v.abs() > 1e-12 {
            let mut w = format!("{s} = 0");
            if !vars.is_empty() {
                w.push_str(", ");
                w.push_str(&sample::describe(vars, &p));
            }
            w.push_str(&format!(" gives rate = {}", sample::fmt_num(v)));
            return Some(w);
        }
    }
    None
}

/// RL009 and RL010.
fn declared(ctx: &mut Ctx, r: &Reaction) {
    let view = ctx.view;
    let Some(math) = r.kinetic_law.as_ref().and_then(|k| k.math.as_ref()) else {
        return;
    };
    let locals = locals(r);
    let direct = idents(&view.inline_functions(math));
    let listed: HashSet<&str> = r
        .reactants
        .iter()
        .chain(&r.products)
        .map(|s| s.species.as_str())
        .chain(r.modifiers.iter().map(|m| m.species.as_str()))
        .collect();
    let mut undeclared: Vec<&String> = direct
        .iter()
        .filter(|id| {
            !locals.contains_key(id.as_str())
                && view.kinds.get(id.as_str()) == Some(&Kind::Species)
                && !listed.contains(id.as_str())
        })
        .collect();
    undeclared.sort();
    for s in undeclared {
        ctx.add(
            RL_UNDECLARED_SPECIES,
            Severity::Warning,
            Element::Reaction(r.id.clone()),
            format!("The rate law uses species {s}, which is not a reactant, product or modifier of the reaction. Declare it as a modifier."),
            None,
        );
    }
    let full: HashSet<String> = expanded_rate(view, r)
        .map(|m| free_symbols(&m).into_iter().collect())
        .unwrap_or_default();
    for m in &r.modifiers {
        if !full.contains(&m.species) {
            ctx.add(
                RL_UNUSED_MODIFIER,
                Severity::Info,
                Element::Reaction(r.id.clone()),
                format!("{} is declared as a modifier but the rate law does not depend on it.", m.species),
                None,
            );
        }
    }
}

fn list(items: &[&String]) -> String {
    match items {
        [] => String::new(),
        [a] => a.to_string(),
        [a, b] => format!("{a} and {b}"),
        _ => {
            let (last, rest) = items.split_last().unwrap();
            format!(
                "{} and {last}",
                rest.iter().map(|s| s.as_str()).collect::<Vec<_>>().join(", ")
            )
        }
    }
}
