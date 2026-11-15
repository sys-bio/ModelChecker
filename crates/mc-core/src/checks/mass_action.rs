//! RL013: a mass-action rate law whose order in a species differs from the
//! species' stoichiometry.
//!
//! For an elementary reaction `2 A + B -> C`, mass action gives
//! `k*A^2*B`: each reactant to the power of its stoichiometry, and the
//! same for products in a reverse term. A rate law of that shape (one
//! forward term, at most one reverse term, each a constant times powers of
//! species) whose powers differ, such as `k*A*B` here, is often a slip.
//! It may be intended for a non-elementary reaction, so this is a warning.
//!
//! Rate laws whose SBO term says mass action are checked more fully under
//! SO007 and skipped here.

use std::collections::BTreeMap;

use sbml_io::{Math, Op};

use super::Ctx;
use super::poly::{Mono, monomials};
use super::scope::{Scope, free_symbols};
use crate::codes::RL_MASS_ACTION_ORDER;
use crate::finding::{Element, Severity};
use crate::fmt::infix;
use crate::sample::fmt_num;
use crate::sbo;
use crate::view::Kind;

pub fn check(ctx: &mut Ctx) {
    let view = ctx.view;
    for r in &view.model.reactions {
        let Some(k) = &r.kinetic_law else { continue };
        let Some(math) = &k.math else { continue };
        let sbo_mass_action = k
            .meta
            .sbo_term
            .as_deref()
            .and_then(sbo::parse)
            .is_some_and(|t| sbo::is_a(t, sbo::MASS_ACTION));
        if sbo_mass_action {
            continue;
        }
        let written = view.inline_functions(math);
        let is_species_id = |id: &str| view.kinds.get(id) == Some(&Kind::Species);
        // Power-law (GMA, S-system) kinetics: species raised to fitted
        // orders, which are not meant to be stoichiometries.
        if power_law(&written, &is_species_id) {
            continue;
        }
        let mut scope = Scope::new(view);
        scope.locals = super::locals(r);
        let value = |m: &sbml_io::Math| scope.at(m, &[], &[]).0;
        let is_species = |id: &str| view.kinds.get(id) == Some(&Kind::Species);
        let Some(monos) = monomials(&written, &is_species, &value) else {
            continue;
        };
        let fwd: Vec<&Mono> = monos.iter().filter(|m| m.0 > 0.0).collect();
        let rev: Vec<&Mono> = monos.iter().filter(|m| m.0 < 0.0).collect();
        if fwd.len() != 1 || rev.len() > 1 {
            continue;
        }
        let mut problems: Vec<String> = Vec::new();
        let reactants = stoichiometry(&r.reactants);
        let products = stoichiometry(&r.products);
        // Only plausibly elementary reactions: whole-number stoichiometries
        // (fractions usually encode volume ratios) and at most three
        // molecules a side (more is a lumped reaction).
        let elementary = |side: &BTreeMap<&str, f64>| {
            side.values().all(|s| s.fract() == 0.0) && side.values().sum::<f64>() <= 3.0
        };
        if !elementary(&reactants) || !elementary(&products) {
            continue;
        }
        // The forward term must be a reactant term; otherwise this is a
        // synthesis or growth law (`mu*m*(1 - m/K)`), not mass action.
        if !reactants.keys().any(|s| fwd[0].1.contains_key(*s)) {
            continue;
        }
        compare(&mut problems, fwd[0], &reactants, "reactant");
        if let [b] = rev.as_slice() {
            compare(&mut problems, b, &products, "product");
        }
        if problems.is_empty() {
            continue;
        }
        ctx.add(
            RL_MASS_ACTION_ORDER,
            Severity::Warning,
            Element::Reaction(r.id.clone()),
            format!(
                "The rate law has the form of mass action, but its powers differ from the \
                 stoichiometry: {}. For an elementary reaction each species appears to the \
                 power of its stoichiometry; if the reaction is not elementary this may be \
                 intended.\nRate law: {}",
                problems.join("; "),
                infix(&written)
            ),
            None,
        );
    }
}

/// Whether `m` raises a species to a power that is not a whole number
/// written as a number: `S^n`, `S^0.5`.
fn power_law(m: &Math, is_species: &dyn Fn(&str) -> bool) -> bool {
    match m {
        Math::Apply(Op::Power, args) if args.len() == 2 => {
            let base_species = free_symbols(&args[0]).iter().any(|s| is_species(s));
            let whole = matches!(&args[1], Math::Number { value, .. } if value.value().fract() == 0.0);
            (base_species && !whole) || args.iter().any(|a| power_law(a, is_species))
        }
        Math::Root { arg, .. } => {
            free_symbols(arg).iter().any(|s| is_species(s)) || power_law(arg, is_species)
        }
        Math::Apply(_, args) | Math::Call(_, args) | Math::CallSymbol(_, args) => {
            args.iter().any(|a| power_law(a, is_species))
        }
        Math::Semantics { math, .. } => power_law(math, is_species),
        _ => false,
    }
}

fn stoichiometry(refs: &[sbml_io::SpeciesReference]) -> BTreeMap<&str, f64> {
    let mut m = BTreeMap::new();
    for sr in refs {
        *m.entry(sr.species.as_str()).or_insert(0.0) += sr.stoichiometry.unwrap_or(1.0);
    }
    m
}

/// Each participant that appears in the term should do so to the power of
/// its stoichiometry. A participant missing from the term altogether is
/// RL006's business (a reactant) or simply irreversibility (a product).
fn compare(problems: &mut Vec<String>, term: &Mono, participants: &BTreeMap<&str, f64>, role: &str) {
    // A term with none of this side's species is not that side's term
    // (a zero-order forward rate, a reverse term on modifiers).
    if !participants.keys().any(|s| term.1.contains_key(*s)) {
        return;
    }
    for (s, st) in participants {
        let Some(e) = term.1.get(*s) else { continue };
        if (e - st).abs() > 1e-9 {
            problems.push(format!(
                "{role} {s} has stoichiometry {} but appears to the power {}",
                fmt_num(*st),
                fmt_num(*e)
            ));
        }
    }
}
