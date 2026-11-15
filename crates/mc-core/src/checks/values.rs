//! Values: parameters in roles they cannot play, reactions switched off,
//! implausible starting values, a wide spread of rate constants, duplicate
//! reactions, and conditions that are always or never true.
//!
//! - VA001 a parameter whose value does not fit its role in a rate law:
//!   a negative rate constant, a saturation or binding constant <= 0, a
//!   Hill coefficient <= 0. The role is read from where the parameter sits:
//!   a factor of the rate is a rate constant (or maximum rate); a constant
//!   added to a concentration in a denominator (`Km + S`) a saturation
//!   constant; a constant dividing a concentration (`S/Ki`) a binding
//!   constant; an exponent on a concentration (`S^n`) a Hill coefficient.
//! - VA002 a rate that is zero for every state: the reaction never runs.
//! - VA003 negative starting amounts, compartments of size <= 0, and values
//!   so large or small they suggest a units slip.
//! - VA004 rate constants spanning ten or more orders of magnitude.
//! - VA005 reactions identical to another (participants and rate law).
//! - VA006 a piecewise condition on the model's state that is always true
//!   or never true, because species cannot be negative.

use std::collections::{BTreeMap, HashSet};

use sbml_io::{Math, Op};

use super::scope::{Scope, free_symbols};
use super::{Ctx, expanded_rate, locals};
use crate::codes::*;
use crate::finding::{Element, Severity};
use crate::fmt::infix;
use crate::sample::fmt_num;
use crate::value::{Tri, Value};
use crate::view::Kind;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Role {
    RateConstant,
    Saturation,
    Binding,
    Hill,
}

pub fn check(ctx: &mut Ctx) {
    roles_and_spread(ctx);
    switched_off(ctx);
    starting_values(ctx);
    duplicates(ctx);
    conditions(ctx);
}

// --------------------------------------------------------- roles

/// Whether `m` depends on a species.
fn has_species(m: &Math, is_species: &dyn Fn(&str) -> bool) -> bool {
    free_symbols(m).iter().any(|s| is_species(s))
}

fn ident(m: &Math) -> Option<&str> {
    match m {
        Math::Ident(id) => Some(id),
        Math::Semantics { math, .. } => ident(math),
        _ => None,
    }
}

/// The factors of a product, through nested products and numerators.
fn factors<'a>(m: &'a Math, out: &mut Vec<&'a Math>) {
    match m {
        Math::Apply(Op::Times, args) => args.iter().for_each(|a| factors(a, out)),
        Math::Apply(Op::Divide, args) if args.len() == 2 => factors(&args[0], out),
        Math::Apply(Op::Minus, args) if args.len() == 1 => factors(&args[0], out),
        Math::Semantics { math, .. } => factors(math, out),
        _ => out.push(m),
    }
}

/// The terms of a sum, through nested sums and differences.
fn terms<'a>(m: &'a Math, out: &mut Vec<&'a Math>) {
    match m {
        Math::Apply(Op::Plus, args) => args.iter().for_each(|a| terms(a, out)),
        Math::Apply(Op::Minus, args) => args.iter().for_each(|a| terms(a, out)),
        Math::Semantics { math, .. } => terms(math, out),
        _ => out.push(m),
    }
}

/// Roles found anywhere in `m` (`den`: inside a denominator).
fn scan(m: &Math, den: bool, is_species: &dyn Fn(&str) -> bool, out: &mut Vec<(String, Role)>) {
    match m {
        Math::Apply(Op::Divide, args) if args.len() == 2 => {
            if let Some(p) = ident(&args[1])
                && has_species(&args[0], is_species)
            {
                out.push((p.to_string(), Role::Binding));
            }
            scan(&args[0], den, is_species, out);
            scan(&args[1], true, is_species, out);
        }
        Math::Apply(Op::Plus, args) => {
            if den && args.iter().any(|a| has_species(a, is_species)) {
                for a in args {
                    if let Some(p) = ident(a) {
                        out.push((p.to_string(), Role::Saturation));
                    }
                    if let Math::Apply(Op::Power, pa) = a
                        && pa.len() == 2
                        && let Some(p) = ident(&pa[0])
                    {
                        out.push((p.to_string(), Role::Saturation));
                    }
                }
            }
            args.iter().for_each(|a| scan(a, den, is_species, out));
        }
        Math::Apply(Op::Power, args) if args.len() == 2 => {
            // A Hill coefficient sits in a saturating denominator
            // (K^n + S^n). An exponent outside one is a kinetic order of a
            // power-law model, which may well be negative (inhibition).
            if den
                && let Some(n) = ident(&args[1])
                && has_species(&args[0], is_species)
            {
                out.push((n.to_string(), Role::Hill));
            }
            args.iter().for_each(|a| scan(a, den, is_species, out));
        }
        Math::Apply(_, args) | Math::Call(_, args) | Math::CallSymbol(_, args) => {
            args.iter().for_each(|a| scan(a, den, is_species, out))
        }
        Math::Piecewise { pieces, otherwise } => {
            for (v, _) in pieces {
                scan(v, den, is_species, out);
            }
            if let Some(o) = otherwise {
                scan(o, den, is_species, out);
            }
        }
        Math::Root { arg, .. } | Math::Log { arg, .. } => scan(arg, den, is_species, out),
        Math::Semantics { math, .. } => scan(math, den, is_species, out),
        _ => {}
    }
}

/// VA001 and VA004.
fn roles_and_spread(ctx: &mut Ctx) {
    let view = ctx.view;
    let is_species = |id: &str| view.kinds.get(id) == Some(&Kind::Species);
    // Rate constants with their values, for the spread.
    let mut constants: Vec<(String, f64)> = Vec::new();
    for r in &view.model.reactions {
        let Some(math) = r.kinetic_law.as_ref().and_then(|k| k.math.as_ref()) else { continue };
        let law = view.inline_functions(math);
        let loc = locals(r);
        let mut found: Vec<(String, Role)> = Vec::new();
        // Factors of each term: rate constants and maximum rates.
        let mut ts = Vec::new();
        terms(&law, &mut ts);
        for t in ts {
            let mut fs = Vec::new();
            factors(t, &mut fs);
            for f in fs {
                if let Some(p) = ident(f) {
                    found.push((p.to_string(), Role::RateConstant));
                }
            }
        }
        scan(&law, false, &is_species, &mut found);
        found.sort();
        found.dedup();
        // A parameter in several roles is judged by the strictest.
        let mut by_param: BTreeMap<String, Role> = BTreeMap::new();
        for (p, role) in found {
            let e = by_param.entry(p).or_insert(role);
            if role != Role::RateConstant {
                *e = role;
            }
        }
        for (p, role) in by_param {
            let value = match loc.get(&p) {
                Some(v) => *v,
                None => {
                    // Only constants: a value a rule or initial assignment
                    // works out (an ionic current, say) has its own meaning.
                    if view.kinds.get(p.as_str()) != Some(&Kind::Parameter)
                        || view.assignment_rules.contains_key(p.as_str())
                        || view.initial_assignments.contains_key(p.as_str())
                    {
                        continue;
                    }
                    match view.range(&p) {
                        Some(rg) if !rg.varies => rg.nominal,
                        _ => continue,
                    }
                }
            };
            if !value.is_finite() {
                continue;
            }
            if role == Role::RateConstant && value > 0.0 {
                constants.push((p.clone(), value));
            }
            let (bad, what, should) = match role {
                Role::RateConstant => (value < 0.0, "multiplies the rate (a rate constant or maximum rate)", "should not be negative"),
                Role::Saturation => (value < 0.0, "is added to a concentration in a denominator (a half-saturation or binding constant)", "should not be negative"),
                Role::Binding => (value <= 0.0, "divides a concentration (a binding or inhibition constant)", "should be positive"),
                Role::Hill => (value <= 0.0, "is the power a concentration is raised to (a Hill coefficient)", "should be positive"),
            };
            if bad {
                let title = match role {
                    Role::RateConstant => format!("Negative rate constant {p}"),
                    Role::Saturation | Role::Binding => format!("Constant {p} is not positive"),
                    Role::Hill => format!("Hill coefficient {p} is not positive"),
                };
                ctx.add_titled(
                    VA_ROLE,
                    Severity::Warning,
                    Element::Reaction(r.id.clone()),
                    title,
                    format!("In this rate law `{p}` {what}, but it is {}; it {should}.", fmt_num(value)),
                    None,
                );
            }
        }
    }

    // VA004: the spread of rate constants.
    constants.sort_by(|a, b| a.0.cmp(&b.0));
    constants.dedup_by(|a, b| a.0 == b.0);
    let (Some(lo), Some(hi)) = (
        constants.iter().min_by(|a, b| a.1.total_cmp(&b.1)),
        constants.iter().max_by(|a, b| a.1.total_cmp(&b.1)),
    ) else {
        return;
    };
    let orders = (hi.1 / lo.1).log10();
    if orders >= 10.0 {
        ctx.add(
            VA_SPREAD,
            Severity::Info,
            Element::Model,
            format!(
                "The rate constants span {:.0} orders of magnitude, from {} = {} to {} = {}. \
                 Processes on very different time scales make a model stiff: it needs a stiff \
                 solver, and some simulators will be slow or fail. If the spread is not real, \
                 a unit may have slipped.",
                orders,
                lo.0,
                fmt_num(lo.1),
                hi.0,
                fmt_num(hi.1)
            ),
            None,
        );
    }
}

// ------------------------------------------------------- switched off

/// VA002. Reactions switched off on purpose (a rate constant set to zero
/// in a model variant) are common, so this is information, and many are
/// listed together.
fn switched_off(ctx: &mut Ctx) {
    let view = ctx.view;
    let mut off: Vec<(String, String, String)> = Vec::new(); // (reaction, why, law)
    for r in &view.model.reactions {
        let Some(law) = expanded_rate(view, r) else { continue };
        let mut scope = Scope::new(view);
        scope.locals = locals(r);
        let (iv, _) = scope.interval(&law);
        if !(iv.lo == 0.0 && iv.hi == 0.0) {
            continue;
        }
        // What is zero: fixed symbols in the rate law with value 0.
        let mut zeros: Vec<String> = free_symbols(&law)
            .into_iter()
            .filter(|s| scope.range(s).is_some_and(|rg| !rg.varies && rg.nominal == 0.0))
            .collect();
        zeros.sort();
        zeros.dedup();
        let why = if zeros.is_empty() {
            String::new()
        } else {
            format!(" because {} {} zero", zeros.join(", "), if zeros.len() == 1 { "is" } else { "are" })
        };
        off.push((r.id.clone(), why, infix(&law)));
    }
    if off.len() >= 5 {
        let list: Vec<String> = off
            .iter()
            .take(15)
            .map(|(r, why, _)| format!("{r}{}", why.replacen(" because", ":", 1)))
            .collect();
        ctx.add(
            VA_SWITCHED_OFF,
            Severity::Info,
            Element::Model,
            format!(
                "{} reactions have a rate of zero whatever the state of the model, so they \
                 never run: {}{}. This is often deliberate (a process switched off in one \
                 version of a model); if not, a value is probably missing.",
                off.len(),
                list.join("; "),
                if off.len() > 15 { format!("; and {} more", off.len() - 15) } else { String::new() }
            ),
            None,
        );
        return;
    }
    for (r, why, law) in off {
        ctx.add(
            VA_SWITCHED_OFF,
            Severity::Info,
            Element::Reaction(r),
            format!(
                "The rate is zero whatever the state of the model{why}, so the reaction never \
                 runs. This may be deliberate (a process switched off); if not, a value is \
                 probably missing.\nRate law: {law}"
            ),
            None,
        );
    }
}

// ----------------------------------------------------- starting values

/// VA003.
fn starting_values(ctx: &mut Ctx) {
    let m = ctx.view.model;
    for s in &m.species {
        // One finding per species, even with both an amount and a
        // concentration. A tiny negative (-2e-27) is a solver's rounding,
        // saved as a starting value: zero in effect, and harmless.
        let negative = [(s.initial_concentration, "concentration"), (s.initial_amount, "amount")]
            .into_iter()
            .find_map(|(v, what)| v.filter(|v| *v < -1e-12).map(|v| (v, what)));
        if let Some((v, what)) = negative {
            ctx.add(
                VA_START,
                Severity::Warning,
                Element::Species(s.id.clone()),
                format!(
                    "{} starts at a negative {what}, {}. An amount of substance cannot be \
                     negative, unless this species stands for something else (a voltage, a \
                     deviation from a reference).",
                    s.id,
                    fmt_num(v)
                ),
                None,
            );
        }
    }
    for c in &m.compartments {
        if let Some(v) = c.size
            && v <= 0.0
            && c.spatial_dimensions != Some(0.0)
        {
            ctx.add(
                VA_START,
                Severity::Error,
                Element::Compartment(c.id.clone()),
                format!(
                    "The compartment's size is {}. A size must be positive: the concentrations \
                     of the species in it cannot be worked out.",
                    fmt_num(v)
                ),
                None,
            );
        }
    }
    // Values so far from ordinary that a unit has probably slipped.
    let extreme = |v: f64| v != 0.0 && v.is_finite() && (v.abs() > 1e25 || v.abs() < 1e-25);
    let mut odd: Vec<String> = Vec::new();
    for s in &m.species {
        for v in [s.initial_concentration, s.initial_amount].into_iter().flatten() {
            if extreme(v) {
                odd.push(format!("{} = {}", s.id, fmt_num(v)));
            }
        }
    }
    for p in &m.parameters {
        if let Some(v) = p.value
            && extreme(v)
        {
            odd.push(format!("{} = {}", p.id, fmt_num(v)));
        }
    }
    for c in &m.compartments {
        if let Some(v) = c.size
            && extreme(v)
        {
            odd.push(format!("{} = {}", c.id, fmt_num(v)));
        }
    }
    if !odd.is_empty() {
        let n = odd.len();
        let shown = odd.iter().take(10).cloned().collect::<Vec<_>>().join(", ");
        ctx.add(
            VA_START,
            Severity::Info,
            Element::Model,
            format!(
                "{} {} extremely large or small (beyond 1e25 or below 1e-25): {shown}{}. Such \
                 values are sometimes right (Avogadro's number, counts of molecules) but often \
                 mean a unit has slipped.",
                n,
                if n == 1 { "value is" } else { "values are" },
                if n > 10 { format!(" and {} more", n - 10) } else { String::new() }
            ),
            None,
        );
    }
}

// ----------------------------------------------------------- duplicates

/// VA005.
fn duplicates(ctx: &mut Ctx) {
    let view = ctx.view;
    let side = |refs: &[sbml_io::SpeciesReference]| -> Vec<(String, String)> {
        let mut v: Vec<(String, String)> = refs
            .iter()
            .map(|r| (r.species.clone(), fmt_num(r.stoichiometry.unwrap_or(1.0))))
            .collect();
        v.sort();
        v
    };
    // (reactants, products, rate law) to the first reaction with them.
    type Side = Vec<(String, String)>;
    let mut seen: BTreeMap<(Side, Side, String), String> = BTreeMap::new();
    for r in &view.model.reactions {
        let Some(law) = expanded_rate(view, r) else { continue };
        // Local parameters stand for their values, so two reactions with
        // their own `k` of the same value are the same.
        let text = infix(&with_values(&law, &locals(r)));
        let key = (side(&r.reactants), side(&r.products), text);
        if key.0.is_empty() && key.1.is_empty() {
            continue;
        }
        match seen.get(&key) {
            Some(first) => ctx.add(
                VA_DUPLICATE,
                Severity::Warning,
                Element::Reaction(r.id.clone()),
                format!(
                    "The reaction is the same as {first}: the same reactants, products and rate \
                     law, so the process is counted twice. Usually a copy-and-paste slip."
                ),
                None,
            ),
            None => {
                seen.insert(key, r.id.clone());
            }
        }
    }
}

/// `m` with the named symbols replaced by their values.
fn with_values(m: &Math, values: &std::collections::HashMap<String, f64>) -> Math {
    let rec = |x: &Math| with_values(x, values);
    match m {
        Math::Ident(id) => match values.get(id) {
            Some(v) if v.is_finite() => Math::real(*v),
            _ => m.clone(),
        },
        Math::Apply(op, a) => Math::Apply(*op, a.iter().map(rec).collect()),
        Math::Call(f, a) => Math::Call(f.clone(), a.iter().map(rec).collect()),
        Math::CallSymbol(s, a) => Math::CallSymbol(s.clone(), a.iter().map(rec).collect()),
        Math::Root { degree, arg } => Math::Root {
            degree: degree.as_ref().map(|d| Box::new(rec(d))),
            arg: Box::new(rec(arg)),
        },
        Math::Log { base, arg } => Math::Log {
            base: base.as_ref().map(|d| Box::new(rec(d))),
            arg: Box::new(rec(arg)),
        },
        Math::Piecewise { pieces, otherwise } => Math::Piecewise {
            pieces: pieces.iter().map(|(v, c)| (rec(v), rec(c))).collect(),
            otherwise: otherwise.as_ref().map(|o| Box::new(rec(o))),
        },
        Math::Semantics { math, .. } => rec(math),
        _ => m.clone(),
    }
}

// ------------------------------------------------------------ conditions

/// The piecewise conditions in `m`, each with whether it is the last piece
/// (nothing after it).
fn piecewise_conditions<'a>(m: &'a Math, out: &mut Vec<(&'a Math, bool)>) {
    match m {
        Math::Piecewise { pieces, otherwise } => {
            for (i, (v, c)) in pieces.iter().enumerate() {
                let last = i + 1 == pieces.len() && otherwise.is_none();
                out.push((c, last));
                piecewise_conditions(v, out);
                piecewise_conditions(c, out);
            }
            if let Some(o) = otherwise {
                piecewise_conditions(o, out);
            }
        }
        Math::Apply(_, a) | Math::Call(_, a) | Math::CallSymbol(_, a) => {
            a.iter().for_each(|x| piecewise_conditions(x, out))
        }
        Math::Root { arg, .. } | Math::Log { arg, .. } => piecewise_conditions(arg, out),
        Math::Semantics { math, .. } => piecewise_conditions(math, out),
        Math::Lambda { body, .. } => piecewise_conditions(body, out),
        _ => {}
    }
}

/// VA006.
fn conditions(ctx: &mut Ctx) {
    let view = ctx.view;
    let mut exprs: Vec<(Element, Math, std::collections::HashMap<String, f64>)> = Vec::new();
    for r in &view.model.reactions {
        if let Some(law) = expanded_rate(view, r) {
            exprs.push((Element::Reaction(r.id.clone()), law, locals(r)));
        }
    }
    for rule in &view.model.rules {
        if let (Some(v), Some(m)) = (rule.variable(), rule.math()) {
            exprs.push((Element::Rule(v.to_string()), view.expand(m, &HashSet::new()), Default::default()));
        }
    }
    for (el, m, loc) in exprs {
        let mut scope = Scope::new(view);
        scope.locals = loc;
        let mut conds = Vec::new();
        piecewise_conditions(&m, &mut conds);
        let mut reported: HashSet<String> = HashSet::new();
        for (c, last) in conds {
            // Only conditions on the model's state: a condition on fixed
            // values is a deliberate switch between model variants.
            if scope.vars(c).is_empty() {
                continue;
            }
            let (v, _) = scope.interval(c);
            // Report only what non-negativity decides: if the condition is
            // decided even with every varying value free to be negative, the
            // model's fixed values decide it, a deliberate setting.
            // (Time stays non-negative: a condition only time decides is
            // decided by the fixed times it is compared with.)
            let free = |id: &str| -> Option<crate::value::Interval> {
                let rg = scope.range(id)?;
                Some(if rg.varies && id != crate::eval::TIME {
                    crate::value::Interval::ALL
                } else {
                    rg.interval
                })
            };
            let anywhere = crate::eval::Eval::new(&free).eval(c);
            if anywhere.truth() != Tri::Unknown {
                continue;
            }
            let text = infix(c);
            let message = match v.truth() {
                Tri::False => format!(
                    "The condition `{text}` can never be true, because species are never \
                     negative, so the branch that depends on it is never used. If it guards \
                     against a case that cannot happen, this is harmless."
                ),
                Tri::True if !last => format!(
                    "The condition `{text}` is always true, because species are never \
                     negative, so the branches after it are never used."
                ),
                _ => continue,
            };
            if reported.insert(text) {
                ctx.add(VA_CONDITION, Severity::Info, el.clone(), message, None);
            }
        }
    }
}
