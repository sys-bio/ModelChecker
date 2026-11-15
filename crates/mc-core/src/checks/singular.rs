//! RL001, RL002: singularities (division by zero, log of a non-positive
//! number, ...) at states a simulation can reach.

use std::collections::HashSet;

use sbml_io::{Math, Rule};

use super::scope::{Scope, search_hazards};
use super::{Ctx, expanded_rate, locals};
use crate::codes::{RL_POSSIBLY_SINGULAR, RL_SINGULAR, RL_SINGULAR_INITIAL, SY_UNSUPPORTED_MATH};
use crate::finding::{Element, Severity};
use crate::fmt::infix;
use crate::value::{HazardKind, Status};
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
        let found = examine(
            ctx,
            &scope,
            &expanded,
            Element::Rule(var.to_string()),
            what,
            &HashSet::new(),
        );
        from_rules.extend(found);
    }

    for r in &view.model.reactions {
        let Some(expanded) = expanded_rate(view, r) else {
            continue;
        };
        let mut scope = Scope::new(view);
        scope.locals = locals(r);
        examine(
            ctx,
            &scope,
            &expanded,
            Element::Reaction(r.id.clone()),
            "rate law",
            &from_rules,
        );
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
        let problem = plain(w.hazard.kind, what, &w.hazard.operand);
        let (code, severity, title, message) = if w.at_initial {
            (
                RL_SINGULAR_INITIAL,
                Severity::Error,
                format!("{} at the start", short(w.hazard.kind)),
                format!(
                    "At the model's starting values, {problem}. The rate cannot be calculated, \
                     so a simulation will fail immediately."
                ),
            )
        } else {
            (
                RL_SINGULAR,
                Severity::Warning,
                format!("{} possible", short(w.hazard.kind)),
                format!(
                    "If the model reaches values like those below, {problem}, and the {what} \
                     cannot be calculated there. Whether a simulation actually gets there depends \
                     on how the model behaves; some models (power-law models, for example) are \
                     written this way on purpose."
                ),
            )
        };
        ctx.add_titled(
            code,
            severity,
            element.clone(),
            title,
            format!("{message}\nExpression: {}", infix(m)),
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
        ctx.add_titled(
            RL_POSSIBLY_SINGULAR,
            Severity::Info,
            element.clone(),
            format!("{} cannot be ruled out", short(h.kind)),
            format!(
                "The checker could not rule out that {}, for some values of the species. It \
                 found no actual values where this happens, so it may well be fine.\nExpression: {}",
                plain(h.kind, what, &h.operand),
                infix(m)
            ),
            None,
        );
    }
    found
}

/// A short name for the problem, for titles.
fn short(kind: HazardKind) -> &'static str {
    match kind {
        // x^-1 is 1/x: to a biologist, dividing by zero.
        HazardKind::DivideByZero | HazardKind::ZeroToNegativePower => "Division by zero",
        HazardKind::LogOfNonPositive => "Logarithm of zero or a negative number",
        HazardKind::RootOfNegative => "Square root of a negative number",
        HazardKind::NegativeBaseFractionalPower => "Negative number to a fractional power",
        HazardKind::Domain => "Function given an impossible value",
    }
}

/// The problem in plain words, as part of a sentence.
fn plain(kind: HazardKind, what: &str, x: &str) -> String {
    match kind {
        HazardKind::DivideByZero => format!("the {what} divides by `{x}`, which is zero"),
        HazardKind::ZeroToNegativePower => format!(
            "the {what} raises `{x}` to a negative power, which means dividing by it, and `{x}` \
             is zero"
        ),
        HazardKind::LogOfNonPositive => format!(
            "the {what} takes the logarithm of `{x}`, which is zero or negative (logarithms \
             exist only for positive numbers)"
        ),
        HazardKind::RootOfNegative => format!(
            "the {what} takes the square root (or another even root) of `{x}`, which is negative"
        ),
        HazardKind::NegativeBaseFractionalPower => format!(
            "the {what} raises `{x}`, which is negative, to a fractional power such as 0.5, \
             which has no real answer"
        ),
        HazardKind::Domain => format!(
            "the {what} applies a function to `{x}` outside the values it accepts (such as \
             arcsin of a number above 1)"
        ),
    }
}
