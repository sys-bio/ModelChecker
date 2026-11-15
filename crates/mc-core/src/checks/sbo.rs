//! SBO term checks.
//!
//! On every element: SO001 an unknown term, SO002 an obsolete one, SO003 a
//! term from the wrong branch for the element (as libSBML checks).
//!
//! On kinetic laws, the term against the formula: SO004 an "irreversible"
//! term on a rate law that can be negative; SO005 a "reversible" term on
//! one that never can; SO006 a mass-action term on a rate law that is not
//! a sum of products of species powers; SO007 a mass-action term whose
//! orders or reactant and product counts differ from the rate law or the
//! reaction; SO008 an enzymatic or Hill term whose shape or reactant count
//! does not fit.

use std::collections::{BTreeMap, HashSet};

use sbml_io::{Math, Meta, Op, Reaction};

use super::scope::{BUDGET, Scope, free_symbols};
use super::{Ctx, expanded_rate, locals};
use crate::codes::*;
use crate::finding::{Element, Severity};
use crate::fmt::infix;
use crate::sample;
use crate::sbo::{self, Term};
use crate::view::Kind;

pub fn check(ctx: &mut Ctx) {
    let m = ctx.view.model;
    use sbo::*;
    let math = &[MATHEMATICAL_EXPRESSION][..];

    branch(ctx, &m.meta, Element::Model, "the model", &[MODELLING_FRAMEWORK]);
    for f in &m.function_definitions {
        branch(ctx, &f.meta, Element::FunctionDefinition(f.id.clone()), "a function definition", math);
    }
    for c in &m.compartments {
        branch(ctx, &c.meta, Element::Compartment(c.id.clone()), "a compartment", &[PHYSICAL_ENTITY]);
    }
    for s in &m.species {
        // Level 2 Versions 2 and 3 used participant roles for species.
        branch(ctx, &s.meta, Element::Species(s.id.clone()), "a species", &[PHYSICAL_ENTITY, PARTICIPANT_ROLE]);
    }
    for p in &m.parameters {
        branch(ctx, &p.meta, Element::Parameter(p.id.clone()), "a parameter", &[SYSTEMS_DESCRIPTION_PARAMETER]);
    }
    for ia in &m.initial_assignments {
        branch(ctx, &ia.meta, Element::InitialAssignment(ia.symbol.clone()), "an initial assignment", math);
    }
    for (i, r) in m.rules.iter().enumerate() {
        let el = match r.variable() {
            Some(v) => Element::Rule(v.to_string()),
            None => Element::AlgebraicRule(i),
        };
        branch(ctx, r.meta(), el, "a rule", math);
    }
    for c in &m.constraints {
        branch(ctx, &c.meta, Element::Model, "a constraint", math);
    }
    for (i, e) in m.events.iter().enumerate() {
        let el = Element::Event(e.id.clone().unwrap_or_else(|| format!("#{}", i + 1)));
        branch(ctx, &e.meta, el.clone(), "an event", &[OCCURRING_ENTITY]);
        if let Some(t) = &e.trigger {
            branch(ctx, &t.meta, el.clone(), "an event trigger", math);
        }
        if let Some(d) = &e.delay {
            branch(ctx, &d.meta, el.clone(), "an event delay", math);
        }
        for a in &e.assignments {
            branch(ctx, &a.meta, el.clone(), "an event assignment", math);
        }
    }
    for r in &m.reactions {
        let el = Element::Reaction(r.id.clone());
        branch(ctx, &r.meta, el.clone(), "a reaction", &[OCCURRING_ENTITY]);
        for sr in r.reactants.iter().chain(&r.products) {
            branch(ctx, &sr.meta, el.clone(), "a reactant or product", &[PARTICIPANT_ROLE]);
        }
        for md in &r.modifiers {
            branch(ctx, &md.meta, el.clone(), "a modifier", &[MODIFIER]);
        }
        if let Some(k) = &r.kinetic_law {
            for p in &k.parameters {
                let lp = Element::LocalParameter {
                    reaction: r.id.clone(),
                    id: p.id.clone(),
                };
                branch(ctx, &p.meta, lp, "a local parameter", &[SYSTEMS_DESCRIPTION_PARAMETER]);
            }
            if let Some(t) = branch(ctx, &k.meta, el, "a kinetic law", &[RATE_LAW]) {
                rate_law(ctx, r, t);
            }
        }
    }
}

/// SO001-SO003 for one element. Returns the term if it is known and in an
/// allowed branch.
fn branch(ctx: &mut Ctx, meta: &Meta, el: Element, what: &str, allowed: &[u32]) -> Option<Term> {
    let written = meta.sbo_term.as_deref()?;
    let Some(id) = sbo::parse(written) else {
        ctx.add(
            SO_UNKNOWN,
            Severity::Warning,
            el,
            format!("The SBO term `{written}` on {what} is not a valid SBO id."),
            None,
        );
        return None;
    };
    let Some(t) = sbo::term(id) else {
        ctx.add(
            SO_UNKNOWN,
            Severity::Warning,
            el,
            format!("{} on {what} is not a term in the Systems Biology Ontology.", sbo::format(id)),
            None,
        );
        return None;
    };
    if t.obsolete {
        ctx.add(
            SO_OBSOLETE,
            Severity::Info,
            el.clone(),
            format!("{} on {what} is obsolete.", sbo::describe(id)),
            None,
        );
    }
    if !allowed.iter().any(|a| sbo::is_a(id, *a)) {
        let want = allowed
            .iter()
            .map(|a| sbo::describe(*a))
            .collect::<Vec<_>>()
            .join(" or ");
        ctx.add(
            SO_WRONG_BRANCH,
            Severity::Warning,
            el,
            format!(
                "{} on {what} is from the wrong part of the ontology: it should be a kind of {want}.",
                sbo::describe(id)
            ),
            None,
        );
        return None;
    }
    Some(t)
}

// ------------------------------------------------------- rate law shape

/// What a term's name says about the rate law.
#[derive(Debug, Default)]
struct Claims {
    reversible: Option<bool>,
    forward_order: Option<u32>,
    reverse_order: Option<u32>,
    reactants: Option<usize>,
    products: Option<usize>,
}

fn word_number(w: &str) -> Option<u32> {
    Some(match w {
        "zeroth" => 0,
        "first" => 1,
        "second" => 2,
        "third" => 3,
        "one" | "uni" => 1,
        "two" | "bi" => 2,
        "three" | "tri" => 3,
        _ => return None,
    })
}

/// The word before `marker` in `name`, as a number.
fn number_before(name: &str, marker: &str) -> Option<u32> {
    let i = name.find(marker)?;
    word_number(name[..i].split_whitespace().last()?)
}

fn claims(name: &str) -> Claims {
    let n = name.to_ascii_lowercase();
    let mut c = Claims {
        reversible: if n.contains("irreversible") || n.contains("irrreversible") {
            Some(false)
        } else if n.contains("reversible") {
            Some(true)
        } else {
            None
        },
        ..Default::default()
    };
    if n.contains("order forward") {
        c.forward_order = number_before(&n, " order forward");
        c.reverse_order = number_before(&n, " order reverse");
    } else {
        c.forward_order = number_before(&n, " order ");
    }
    c.reactants = number_before(&n, " reactant").map(|x| x as usize);
    c.products = number_before(&n, " product").map(|x| x as usize);
    // Enzyme mechanism names: unireactant, Bi-Uni, two substrates, ...
    if n.contains("unireactant") {
        c.reactants = Some(1);
    } else if n.contains("bireactant") || n.contains("two substrates") {
        c.reactants = Some(2);
    } else if n.contains("trireactant") {
        c.reactants = Some(3);
    }
    for (pat, r, p) in [("uni-uni", 1, 1), ("uni-bi", 1, 2), ("bi-uni", 2, 1), ("bi-bi", 2, 2)] {
        if n.contains(pat) {
            c.reactants = Some(r);
            c.products = Some(p);
        }
    }
    c
}

fn rate_law(ctx: &mut Ctx, r: &Reaction, t: Term) {
    let view = ctx.view;
    let Some(law) = expanded_rate(view, r) else { return };
    // The rate law as written, with only function definitions inlined: its
    // shape is what the SBO term describes, and expanding assignment rules
    // would hide the species they stand for.
    let Some(written) = r.kinetic_law.as_ref().and_then(|k| k.math.as_ref()) else {
        return;
    };
    let written = view.inline_functions(written);
    let el = Element::Reaction(r.id.clone());
    let named = sbo::describe(t.id);
    let c = claims(t.name);
    let mut scope = Scope::new(view);
    scope.locals = locals(r);

    // Reversibility.
    match c.reversible {
        Some(false) => {
            if let Some(w) = negative_point(&scope, &law) {
                ctx.add(
                    SO_IRREVERSIBLE_NEGATIVE,
                    Severity::Warning,
                    el.clone(),
                    format!(
                        "The kinetic law is marked {named}, but the rate can be negative, so \
                         the reaction can run backwards.\nRate law: {}",
                        infix(&law)
                    ),
                    Some(w),
                );
            }
        }
        Some(true) => {
            let (iv, _) = scope.interval(&law);
            if iv.lo >= 0.0 {
                ctx.add(
                    SO_REVERSIBLE_NEVER_NEGATIVE,
                    Severity::Warning,
                    el.clone(),
                    format!(
                        "The kinetic law is marked {named}, but the rate can never be \
                         negative: the rate law describes an irreversible reaction.\nRate law: {}",
                        infix(&law)
                    ),
                    None,
                );
            }
        }
        None => {}
    }

    let reactants = participants(&r.reactants);
    let products = participants(&r.products);
    let modifiers: HashSet<&str> = r.modifiers.iter().map(|m| m.species.as_str()).collect();
    let is_species = |id: &str| view.kinds.get(id) == Some(&Kind::Species);

    // Mass action.
    let name = t.name.to_ascii_lowercase();
    let discrete = name.contains("discrete");
    let stimulator = name.contains("stimulator");
    if sbo::is_a(t.id, sbo::MASS_ACTION) && !discrete && !stimulator {
        let value = |m: &Math| {
            let (v, _) = scope.at(m, &[], &[]);
            v
        };
        let Some(monos) = monomials(&written, &is_species, &value) else {
            ctx.add(
                SO_NOT_MASS_ACTION,
                Severity::Warning,
                el,
                format!(
                    "The kinetic law is marked {named}, but the rate law is not mass action: \
                     a sum of terms, each a constant times powers of species.\nRate law: {}",
                    infix(&written)
                ),
                None,
            );
            return;
        };
        let fwd: Vec<&Mono> = monos.iter().filter(|m| m.0 > 0.0).collect();
        let rev: Vec<&Mono> = monos.iter().filter(|m| m.0 < 0.0).collect();
        let mut problems: Vec<String> = Vec::new();
        if fwd.len() != 1 {
            problems.push(format!(
                "it should have one forward (positive) term, but has {}",
                fwd.len()
            ));
        }
        if c.reversible == Some(true) && rev.len() != 1 {
            problems.push(format!(
                "a reversible mass-action rate law should have one reverse (negative) term, but has {}",
                rev.len()
            ));
        }
        if let [f] = fwd.as_slice() {
            let side = Side {
                name: "forward",
                role: "reactant",
                order: c.forward_order,
                count: c.reactants,
                participants: &reactants,
            };
            compare_term(&mut problems, f, &side, &modifiers);
        }
        if let [b] = rev.as_slice() {
            let side = Side {
                name: "reverse",
                role: "product",
                order: c.reverse_order,
                count: c.products,
                participants: &products,
            };
            compare_term(&mut problems, b, &side, &modifiers);
        }
        if !problems.is_empty() {
            ctx.add(
                SO_MASS_ACTION_MISMATCH,
                Severity::Warning,
                el,
                format!(
                    "The kinetic law is marked {named}, but {}.\nRate law: {}",
                    problems.join("; "),
                    infix(&written)
                ),
                None,
            );
        }
        return;
    }

    // Enzymatic and Hill.
    if sbo::is_a(t.id, sbo::MODULAR) {
        return;
    }
    let mut problems: Vec<String> = Vec::new();
    let n_reactants = reactants.len();
    let n_products = products.len();
    if let Some(n) = c.reactants
        && n != n_reactants
    {
        problems.push(format!(
            "the term is for {} reactant{}, but the reaction has {n_reactants}",
            n,
            if n == 1 { "" } else { "s" }
        ));
    }
    if let Some(n) = c.products
        && n != n_products
    {
        problems.push(format!(
            "the term is for {} product{}, but the reaction has {n_products}",
            n,
            if n == 1 { "" } else { "s" }
        ));
    }
    let deps: HashSet<String> = free_symbols(&written).into_iter().collect();
    // A product in the denominator counts too: the substrate is sometimes
    // written as a total minus the product, `(Rt - Rp)/(K + Rt - Rp)`. A
    // modifier alone does not: `v*S*E/(K + E)` saturates in the enzyme.
    if sbo::is_a(t.id, sbo::MICHAELIS_MENTEN_FAMILY)
        && !reactants
            .keys()
            .chain(products.keys())
            .any(|s| saturates(&written, s))
    {
        problems.push(
            "a Michaelis-Menten-type rate law should saturate in its substrate (a substrate \
             in a denominator, as in V*S/(K + S)), but this one does not"
                .to_string(),
        );
    }
    if (sbo::is_a(t.id, sbo::HILL) || sbo::is_a(t.id, sbo::REVERSIBLE_HILL)) && !has_species_power(&written, &is_species) {
        problems.push(
            "a Hill-type rate law should raise a species to a power (S^n), but this one does not"
                .to_string(),
        );
    }
    if sbo::is_a(t.id, sbo::ENZYMATIC_MODULATED) {
        let modulated = deps
            .iter()
            .any(|d| is_species(d) && !reactants.contains_key(d.as_str()));
        if !modulated {
            problems.push(
                "a modulated enzymatic rate law should depend on a modifier (an inhibitor or \
                 activator) or a product, but this one depends only on its reactants"
                    .to_string(),
            );
        }
    }
    if !problems.is_empty() {
        ctx.add(
            SO_SHAPE_MISMATCH,
            Severity::Warning,
            el,
            format!(
                "The kinetic law is marked {named}, but {}.\nRate law: {}",
                problems.join("; "),
                infix(&written)
            ),
            None,
        );
    }
}

/// Species to total stoichiometry.
fn participants(refs: &[sbml_io::SpeciesReference]) -> BTreeMap<&str, f64> {
    let mut m = BTreeMap::new();
    for sr in refs {
        *m.entry(sr.species.as_str()).or_insert(0.0) += sr.stoichiometry.unwrap_or(1.0);
    }
    m
}

/// A reachable point where the rate is negative, as text.
fn negative_point(scope: &Scope, law: &Math) -> Option<String> {
    let (iv, _) = scope.interval(law);
    if iv.lo >= 0.0 {
        return None;
    }
    let vars = scope.vars(law);
    let extras = scope.extras(law);
    for p in sample::points(&vars, &extras, BUDGET) {
        if !scope.reachable(&vars, &p) {
            continue;
        }
        let (v, _) = scope.at(law, &vars, &p);
        if v.is_finite() && v < -1e-12 {
            let mut w = sample::describe(&vars, &p);
            if !w.is_empty() {
                w.push(' ');
            }
            w.push_str(&format!("gives rate = {}", sample::fmt_num(v)));
            return Some(w);
        }
    }
    None
}

/// One side of a mass-action rate law: what the SBO term claims for it,
/// and the reaction's participants there.
struct Side<'a> {
    /// "forward" or "reverse".
    name: &'a str,
    /// "reactant" or "product".
    role: &'a str,
    order: Option<u32>,
    count: Option<usize>,
    participants: &'a BTreeMap<&'a str, f64>,
}

/// Compares one mass-action term with what the SBO term claims and with
/// the reaction's participants on that side.
fn compare_term(problems: &mut Vec<String>, term: &Mono, s: &Side, modifiers: &HashSet<&str>) {
    let (side, role, order, count, participants) = (s.name, s.role, s.order, s.count, s.participants);
    let actual: f64 = term.1.values().sum();
    if let Some(o) = order
        && (actual - o as f64).abs() > 1e-9
    {
        problems.push(format!(
            "the term says {} order {side}, and the {side} part is order {}",
            ordinal(o),
            sample::fmt_num(actual)
        ));
    }
    if let Some(n) = count
        && n != term.1.len()
    {
        problems.push(format!(
            "the term says {n} {role}{}, and the {side} part involves {} species",
            if n == 1 { "" } else { "s" },
            term.1.len()
        ));
    }
    // Against the reaction: each participant to the power of its
    // stoichiometry. Declared modifiers (catalysts, templates) may also
    // appear; the order check above already counts them.
    let expected: BTreeMap<&str, f64> = participants.clone();
    let got: BTreeMap<&str, f64> = term
        .1
        .iter()
        .map(|(k, v)| (k.as_str(), *v))
        .filter(|(k, _)| participants.contains_key(k) || !modifiers.contains(k))
        .collect();
    let same = expected.len() == got.len()
        && expected
            .iter()
            .all(|(k, v)| got.get(k).is_some_and(|g| (g - v).abs() < 1e-9));
    if !same {
        problems.push(format!(
            "the {side} part is {} but the reaction's {role}s call for {}",
            describe_mono(&got),
            describe_mono(&expected)
        ));
    }
}

fn ordinal(o: u32) -> &'static str {
    match o {
        0 => "zeroth",
        1 => "first",
        2 => "second",
        3 => "third",
        _ => "higher",
    }
}

fn describe_mono(m: &BTreeMap<&str, f64>) -> String {
    if m.is_empty() {
        return "a constant".to_string();
    }
    m.iter()
        .map(|(s, e)| {
            if (*e - 1.0).abs() < 1e-9 {
                s.to_string()
            } else {
                format!("{s}^{}", sample::fmt_num(*e))
            }
        })
        .collect::<Vec<_>>()
        .join("*")
}

/// A term of a mass-action expression: its sign and species exponents.
type Mono = (f64, BTreeMap<String, f64>);

const MAX_TERMS: usize = 32;

/// `m` as a sum of signed monomials in the species, if it is one. Factors
/// free of species (rate constants, compartments) are dropped, keeping
/// their sign.
fn monomials(
    m: &Math,
    is_species: &dyn Fn(&str) -> bool,
    value: &dyn Fn(&Math) -> f64,
) -> Option<Vec<Mono>> {
    let species_free = free_symbols(m).iter().all(|s| !is_species(s));
    if species_free {
        let v = value(m);
        let sign = if v < 0.0 { -1.0 } else { 1.0 };
        return Some(vec![(sign, BTreeMap::new())]);
    }
    let rec = |x: &Math| monomials(x, is_species, value);
    let monos = match m {
        Math::Ident(id) => vec![(1.0, BTreeMap::from([(id.clone(), 1.0)]))],
        Math::Semantics { math, .. } => rec(math)?,
        Math::Apply(Op::Plus, args) => {
            let mut out = Vec::new();
            for a in args {
                out.extend(rec(a)?);
            }
            out
        }
        Math::Apply(Op::Minus, args) => match args.as_slice() {
            [a] => negate(rec(a)?),
            [a, b] => {
                let mut out = rec(a)?;
                out.extend(negate(rec(b)?));
                out
            }
            _ => return None,
        },
        Math::Apply(Op::Times, args) => {
            let mut acc: Vec<Mono> = vec![(1.0, BTreeMap::new())];
            for a in args {
                let terms = rec(a)?;
                let mut next = Vec::new();
                for (s1, e1) in &acc {
                    for (s2, e2) in &terms {
                        let mut e = e1.clone();
                        for (k, v) in e2 {
                            *e.entry(k.clone()).or_insert(0.0) += v;
                        }
                        next.push((s1 * s2, e));
                    }
                }
                if next.len() > MAX_TERMS {
                    return None;
                }
                acc = next;
            }
            acc
        }
        Math::Apply(Op::Divide, args) => {
            let [a, b] = args.as_slice() else { return None };
            if free_symbols(b).iter().any(|s| is_species(s)) {
                return None;
            }
            let sign = if value(b) < 0.0 { -1.0 } else { 1.0 };
            rec(a)?.into_iter().map(|(s, e)| (s * sign, e)).collect()
        }
        Math::Apply(Op::Power, args) => {
            let [a, b] = args.as_slice() else { return None };
            if free_symbols(b).iter().any(|s| is_species(s)) {
                return None;
            }
            let p = value(b);
            if !p.is_finite() {
                return None;
            }
            let base = rec(a)?;
            let [(s, e)] = base.as_slice() else { return None };
            let sign = if *s < 0.0 {
                if p.fract() != 0.0 {
                    return None;
                }
                if (p as i64) % 2 == 0 { 1.0 } else { -1.0 }
            } else {
                1.0
            };
            vec![(sign, e.iter().map(|(k, v)| (k.clone(), v * p)).collect())]
        }
        _ => return None,
    };
    // Drop species with exponent 0, then merge like terms of the same sign.
    let mut merged: Vec<Mono> = Vec::new();
    for (s, mut e) in monos {
        e.retain(|_, v| v.abs() > 1e-12);
        if let Some(x) = merged.iter_mut().find(|x| x.1 == e && x.0 == s) {
            x.0 = s;
        } else {
            merged.push((s, e));
        }
    }
    Some(merged)
}

fn negate(v: Vec<Mono>) -> Vec<Mono> {
    v.into_iter().map(|(s, e)| (-s, e)).collect()
}

/// Whether `species` appears in a denominator of `m`.
fn saturates(m: &Math, species: &str) -> bool {
    match m {
        Math::Apply(Op::Divide, args) if args.len() == 2 => {
            free_symbols(&args[1]).iter().any(|s| s == species)
                || args.iter().any(|a| saturates(a, species))
        }
        Math::Apply(Op::Power, args) if args.len() == 2 => {
            // x^-1 is a division too.
            let negative = matches!(&args[1], Math::Number { value, .. } if value.value() < 0.0);
            (negative && free_symbols(&args[0]).iter().any(|s| s == species))
                || args.iter().any(|a| saturates(a, species))
        }
        Math::Apply(_, args) | Math::Call(_, args) | Math::CallSymbol(_, args) => {
            args.iter().any(|a| saturates(a, species))
        }
        Math::Piecewise { pieces, otherwise } => {
            pieces.iter().any(|(v, _)| saturates(v, species))
                || otherwise.as_deref().is_some_and(|o| saturates(o, species))
        }
        Math::Semantics { math, .. } => saturates(math, species),
        _ => false,
    }
}

/// Whether `m` raises something involving a species to a power other
/// than a whole number up to 1 (S^n, S^h, S^2.5).
fn has_species_power(m: &Math, is_species: &dyn Fn(&str) -> bool) -> bool {
    match m {
        Math::Apply(Op::Power, args) if args.len() == 2 => {
            let base_species = free_symbols(&args[0]).iter().any(|s| is_species(s));
            let trivial = matches!(&args[1], Math::Number { value, .. } if value.value() == 1.0);
            (base_species && !trivial) || args.iter().any(|a| has_species_power(a, is_species))
        }
        Math::Apply(_, args) | Math::Call(_, args) | Math::CallSymbol(_, args) => {
            args.iter().any(|a| has_species_power(a, is_species))
        }
        Math::Root { arg, .. } => has_species_power(arg, is_species),
        Math::Piecewise { pieces, otherwise } => {
            pieces.iter().any(|(v, _)| has_species_power(v, is_species))
                || otherwise.as_deref().is_some_and(|o| has_species_power(o, is_species))
        }
        Math::Semantics { math, .. } => has_species_power(math, is_species),
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn claims_from_names() {
        let c = claims("mass action rate law for first order forward, second order reverse, reversible reactions, two products, continuous scheme");
        assert_eq!(c.reversible, Some(true));
        assert_eq!(c.forward_order, Some(1));
        assert_eq!(c.reverse_order, Some(2));
        assert_eq!(c.products, Some(2));
        let c = claims("mass action rate law for second order irreversible reactions, two reactants, continuous scheme");
        assert_eq!(c.reversible, Some(false));
        assert_eq!(c.forward_order, Some(2));
        assert_eq!(c.reactants, Some(2));
        let c = claims("Ordered Bi-Uni mechanism rate law");
        assert_eq!((c.reactants, c.products), (Some(2), Some(1)));
        let c = claims("enzymatic rate law for irreversible non-modulated non-interacting unireactant enzymes");
        assert_eq!((c.reversible, c.reactants), (Some(false), Some(1)));
    }
}
