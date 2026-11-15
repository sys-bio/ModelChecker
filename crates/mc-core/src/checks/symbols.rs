//! Symbol checks: SY001 undefined, SY002-SY004 unused, SY005 shadowing,
//! SY006 no value.

use std::collections::HashSet;

use sbml_io::{Math, Rule};

use super::Ctx;
use crate::codes::*;
use crate::finding::{Element, Severity};

pub fn check(ctx: &mut Ctx) {
    let view = ctx.view;
    let m = view.model;

    // Every expression in the model, with the element it belongs to and the
    // extra names in scope there (local parameters, lambda arguments).
    let mut exprs: Vec<(Element, &Math, HashSet<&str>)> = Vec::new();
    for f in &m.function_definitions {
        if let Some(math) = &f.math {
            let args: HashSet<&str> = match strip(math) {
                Math::Lambda { args, .. } => args.iter().map(String::as_str).collect(),
                _ => HashSet::new(),
            };
            exprs.push((Element::FunctionDefinition(f.id.clone()), math, args));
        }
    }
    for (i, rule) in m.rules.iter().enumerate() {
        if let Some(math) = rule.math() {
            let el = match rule.variable() {
                Some(v) => Element::Rule(v.to_string()),
                None => Element::AlgebraicRule(i),
            };
            exprs.push((el, math, HashSet::new()));
        }
    }
    for ia in &m.initial_assignments {
        if let Some(math) = &ia.math {
            exprs.push((Element::InitialAssignment(ia.symbol.clone()), math, HashSet::new()));
        }
    }
    for (i, e) in m.events.iter().enumerate() {
        let el = Element::Event(e.id.clone().unwrap_or_else(|| format!("#{}", i + 1)));
        let maths = e
            .trigger
            .iter()
            .filter_map(|t| t.math.as_ref())
            .chain(e.delay.iter().filter_map(|d| d.math.as_ref()))
            .chain(e.priority.iter().filter_map(|p| p.math.as_ref()))
            .chain(e.assignments.iter().filter_map(|a| a.math.as_ref()));
        for math in maths {
            exprs.push((el.clone(), math, HashSet::new()));
        }
    }
    for c in &m.constraints {
        if let Some(math) = &c.math {
            exprs.push((Element::Model, math, HashSet::new()));
        }
    }
    for r in &m.reactions {
        let el = Element::Reaction(r.id.clone());
        let locals: HashSet<&str> = r
            .kinetic_law
            .iter()
            .flat_map(|k| &k.parameters)
            .map(|p| p.id.as_str())
            .collect();
        if let Some(math) = r.kinetic_law.as_ref().and_then(|k| k.math.as_ref()) {
            exprs.push((el.clone(), math, locals));
        }
        for sr in r.reactants.iter().chain(&r.products) {
            if let Some(math) = &sr.stoichiometry_math {
                exprs.push((el.clone(), math, HashSet::new()));
            }
        }
    }

    // SY001: undefined. Also gathers what is used, for SY002 and SY004.
    let mut used: HashSet<String> = HashSet::new();
    for (el, math, scope) in &exprs {
        let mut undefined: Vec<String> = Vec::new();
        math.visit_idents(&mut |id| {
            if scope.contains(id) {
                return;
            }
            used.insert(id.to_string());
            if !view.kinds.contains_key(id) && !undefined.iter().any(|u| u == id) {
                undefined.push(id.to_string());
            }
        });
        for id in undefined {
            ctx.add(
                SY_UNDEFINED,
                Severity::Error,
                el.clone(),
                format!("`{id}` is not defined anywhere in the model."),
                None,
            );
        }
    }
    // Ids used as attributes count as used.
    if let Some(cf) = &m.conversion_factor {
        used.insert(cf.clone());
    }
    for s in &m.species {
        if let Some(cf) = &s.conversion_factor {
            used.insert(cf.clone());
        }
    }

    // SY002: global parameters nothing uses and nothing sets.
    let set: HashSet<&str> = m
        .rules
        .iter()
        .filter_map(Rule::variable)
        .chain(view.event_targets.iter().copied())
        .collect();
    for p in &m.parameters {
        if !used.contains(&p.id) && !set.contains(p.id.as_str()) {
            ctx.add(
                SY_UNUSED_PARAMETER,
                Severity::Info,
                Element::Parameter(p.id.clone()),
                "The parameter is never used.",
                None,
            );
        }
    }
    // SY004: functions never called.
    for f in &m.function_definitions {
        if !used.contains(&f.id) {
            ctx.add(
                SY_UNUSED_FUNCTION,
                Severity::Info,
                Element::FunctionDefinition(f.id.clone()),
                "The function is never called.",
                None,
            );
        }
    }

    // SY003, SY005 and local values.
    for r in &m.reactions {
        let Some(k) = &r.kinetic_law else { continue };
        let in_law = k.math.as_ref().map(crate::view::idents).unwrap_or_default();
        for p in &k.parameters {
            let el = Element::LocalParameter {
                reaction: r.id.clone(),
                id: p.id.clone(),
            };
            if !in_law.contains(&p.id) {
                ctx.add(
                    SY_UNUSED_LOCAL,
                    Severity::Info,
                    el.clone(),
                    "The local parameter is not used in the rate law.",
                    None,
                );
            }
            if view.kinds.contains_key(p.id.as_str()) {
                ctx.add(
                    SY_SHADOWED,
                    Severity::Info,
                    el.clone(),
                    format!("The local parameter hides the global {} `{}` inside this rate law.", kind_name(view.kinds[p.id.as_str()]), p.id),
                    None,
                );
            }
            if p.value.is_none() {
                ctx.add(
                    SY_NO_VALUE,
                    Severity::Warning,
                    el,
                    "The local parameter has no value.",
                    None,
                );
            }
        }
    }

    // SY006: global symbols with no value from anywhere.
    let given = |id: &str| {
        view.initial_assignments.contains_key(id)
            || view.assignment_rules.contains_key(id)
            || view.algebraic.contains(id)
    };
    for p in &m.parameters {
        if p.value.is_none() && !given(&p.id) {
            ctx.add(
                SY_NO_VALUE,
                Severity::Warning,
                Element::Parameter(p.id.clone()),
                "The parameter has no value, initial assignment or rule.",
                None,
            );
        }
    }
    for s in &m.species {
        if s.initial_amount.is_none() && s.initial_concentration.is_none() && !given(&s.id) {
            ctx.add(
                SY_NO_VALUE,
                Severity::Warning,
                Element::Species(s.id.clone()),
                "The species has no initial amount or concentration, initial assignment or rule.",
                None,
            );
        }
    }
    // In Level 2 a missing size means 1; in Level 3 it means unknown.
    if view.level >= 3 {
        for c in &m.compartments {
            let dims = c.spatial_dimensions.unwrap_or(3.0);
            if c.size.is_none() && dims != 0.0 && !given(&c.id) {
                ctx.add(
                    SY_NO_VALUE,
                    Severity::Warning,
                    Element::Compartment(c.id.clone()),
                    "The compartment has no size, initial assignment or rule.",
                    None,
                );
            }
        }
    }
}

fn strip(m: &Math) -> &Math {
    match m {
        Math::Semantics { math, .. } => strip(math),
        _ => m,
    }
}

fn kind_name(k: crate::view::Kind) -> &'static str {
    use crate::view::Kind::*;
    match k {
        Species => "species",
        Compartment => "compartment",
        Parameter => "parameter",
        Reaction => "reaction",
        SpeciesReference => "species reference",
        Function => "function",
    }
}
