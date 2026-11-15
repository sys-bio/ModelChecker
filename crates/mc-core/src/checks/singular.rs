//! RL001, RL002: singularities (division by zero, log of a non-positive
//! number, ...) at states a simulation can reach.

use std::collections::HashSet;

use sbml_io::{Math, Rule};

use super::scope::{Scope, search_hazards};
use super::{Ctx, expanded_rate, locals};
use crate::codes::{RL_POSSIBLY_SINGULAR, RL_SINGULAR, RL_SINGULAR_INITIAL, SY_UNSUPPORTED_MATH};
use crate::finding::{Element, Severity};
use crate::fmt::infix;
use crate::value::Status;
use crate::view::has_raw;

pub fn check(ctx: &mut Ctx) {
    let view = ctx.view;
    // Problems found in rules, so an assignment rule's singularity is not
    // reported again in every rate law that uses it.
    let mut from_rules: HashSet<(String, String)> = HashSet::new();

    for rule in &view.model.rules {
        let (Some(var), Some(math)) = (rule.variable(), rule.math()) else {
            continue;
        };
        let what = match rule {
            Rule::Rate { .. } => "rate rule",
            _ => "assignment rule",
        };
        let expanded = view.expand(math, &HashSet::new());
        let scope = Scope::new(view);
        let found = examine(ctx, &scope, &expanded, Element::Rule(var.to_string()), what, &HashSet::new());
        from_rules.extend(found);
    }

    for r in &view.model.reactions {
        let Some(expanded) = expanded_rate(view, r) else {
            continue;
        };
        let mut scope = Scope::new(view);
        scope.locals = locals(r);
        examine(ctx, &scope, &expanded, Element::Reaction(r.id.clone()), "rate law", &from_rules);
    }
}

fn examine(
    ctx: &mut Ctx,
    scope: &Scope,
    m: &Math,
    element: Element,
    what: &str,
    skip: &HashSet<(String, String)>,
) -> HashSet<(String, String)> {
    let mut seen = skip.clone();
    let mut found = HashSet::new();
    if has_raw(m) {
        ctx.add(
            SY_UNSUPPORTED_MATH,
            Severity::Info,
            element.clone(),
            format!("The {what} contains MathML that could not be analysed; results for it may be incomplete."),
            None,
        );
    }
    let key = |kind: crate::value::HazardKind, operand: &str| {
        (kind.describe().to_string(), operand.to_string())
    };

    // Definite problems, with a point that shows them.
    let certain = search_hazards(scope, m);
    for w in &certain {
        let k = key(w.hazard.kind, &w.hazard.operand);
        if !seen.insert(k.clone()) {
            continue;
        }
        found.insert(k);
        let (code, severity, when) = if w.at_initial {
            (RL_SINGULAR_INITIAL, Severity::Error, "at the initial values")
        } else {
            (RL_SINGULAR, Severity::Warning, "at some non-negative state")
        };
        ctx.add(
            code,
            severity,
            element.clone(),
            format!(
                "The {what} has a {} in `{}` {when}.\nExpression: {}",
                w.hazard.kind.describe(),
                w.hazard.operand,
                infix(m)
            ),
            Some(w.witness.clone()),
        );
    }

    // Possible problems no point confirmed. Hazards in a piecewise branch
    // whose condition might be their guard are left out.
    let (_, hazards) = scope.interval(m);
    for h in hazards {
        if h.status == Status::Safe || h.guarded {
            continue;
        }
        if certain
            .iter()
            .any(|w| w.hazard.kind == h.kind && w.hazard.operand == h.operand)
        {
            continue;
        }
        let k = key(h.kind, &h.operand);
        if !seen.insert(k.clone()) {
            continue;
        }
        found.insert(k);
        ctx.add(
            RL_POSSIBLY_SINGULAR,
            Severity::Info,
            element.clone(),
            format!(
                "The {what} may have a {} in `{}` for some non-negative species values.\nExpression: {}",
                h.kind.describe(),
                h.operand,
                infix(m)
            ),
            None,
        );
    }
    found
}
