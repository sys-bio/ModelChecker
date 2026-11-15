//! Network structure.
//!
//! - NW001 only produced, NW002 only consumed, NW003 unused species;
//! - NW004 a network in separate parts that never affect each other;
//! - NW005 conserved totals (conservation laws), from the left null space
//!   of the stoichiometry matrix;
//! - NW006 species that can never become non-zero: they start at zero, and
//!   every reaction that would make one of them needs another of them as a
//!   reactant (an empty siphon). The reactions consuming them never run.
//! - NW007 a reaction that changes no variable species.

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

use sbml_io::Reaction;

use super::Ctx;
use super::scope::free_symbols;
use crate::codes::*;
use crate::finding::{Element, Severity};
use crate::sample::fmt_num;
use crate::view::View;

/// When at least this many species share a finding, they are listed in one
/// model-level finding.
const GROUP_AT: usize = 5;

pub fn check(ctx: &mut Ctx) {
    let view = ctx.view;
    let m = view.model;
    let held = |id: &str| {
        m.species
            .iter()
            .find(|s| s.id == id)
            .is_some_and(|s| s.boundary_condition == Some(true) || s.constant == Some(true))
    };
    // Set other than by reactions: rules, events, algebraic rules.
    let set_otherwise = |id: &str| {
        view.rate_rules.contains(id)
            || view.event_targets.contains(id)
            || view.algebraic.contains(id)
            || view.assignment_rules.contains_key(id)
    };
    let nets: Vec<BTreeMap<String, f64>> = m.reactions.iter().map(net_stoichiometry).collect();
    let used = used_symbols(view);

    // NW001, NW002, NW003.
    let mut only_produced: Vec<(String, bool)> = Vec::new();
    let mut only_consumed: Vec<String> = Vec::new();
    let mut unused: Vec<String> = Vec::new();
    for s in &m.species {
        let id = s.id.as_str();
        if set_otherwise(id) {
            continue;
        }
        let in_reactions = m.reactions.iter().any(|r| {
            r.reactants
                .iter()
                .chain(&r.products)
                .any(|x| x.species == id)
                || r.modifiers.iter().any(|x| x.species == id)
        });
        if !in_reactions && !used.contains(id) {
            unused.push(s.id.clone());
            continue;
        }
        if held(id) {
            continue;
        }
        let produced = nets.iter().any(|n| n.get(id).is_some_and(|x| *x > 0.0));
        let consumed = nets.iter().any(|n| n.get(id).is_some_and(|x| *x < 0.0));
        if produced && !consumed {
            only_produced.push((s.id.clone(), used.contains(id)));
        } else if consumed && !produced {
            only_consumed.push(s.id.clone());
        }
    }
    report_species(
        ctx,
        NW_ONLY_PRODUCED,
        only_produced.iter().map(|(s, _)| s.clone()).collect(),
        |s| {
            let feeds_back = only_produced.iter().any(|(x, u)| x == s && *u);
            format!(
                "{s} is produced but never consumed, so it accumulates. {}",
                if feeds_back {
                    "It is used elsewhere in the model, so its build-up affects the dynamics."
                } else {
                    "Nothing else uses it, so it is an output of the model."
                }
            )
        },
        |list| {
            format!("These species are produced but never consumed, so they accumulate: {list}.")
        },
    );
    report_species(
        ctx,
        NW_ONLY_CONSUMED,
        only_consumed,
        |s| format!("{s} is consumed but never produced, so it runs down over time."),
        |list| {
            format!(
                "These species are consumed but never produced, so they run down over time: {list}."
            )
        },
    );
    report_species(
        ctx,
        NW_UNUSED_SPECIES,
        unused,
        |s| format!("{s} takes part in no reaction, rule or event, and nothing uses it."),
        |list| {
            format!(
                "These species take part in no reaction, rule or event, and nothing uses them: {list}."
            )
        },
    );

    disconnected(ctx);
    conservation(ctx, &nets, &held, &set_otherwise);
    never_present(ctx, &nets, &held, &set_otherwise);

    // NW007.
    let no_effect: Vec<&str> = m
        .reactions
        .iter()
        .zip(&nets)
        .filter(|(r, net)| {
            // A reaction with nothing at all is RL011's.
            !(r.reactants.is_empty() && r.products.is_empty()) && net.keys().all(|s| held(s))
        })
        .map(|(r, _)| r.id.as_str())
        .collect();
    if no_effect.len() >= GROUP_AT {
        ctx.add(
            NW_NO_EFFECT,
            Severity::Info,
            Element::Model,
            format!(
                "{} reactions change no variable species (their species are all boundary or \
                 constant, or appear equally on both sides), so they only carry fluxes: {}.",
                no_effect.len(),
                abbreviate(&no_effect, 12)
            ),
            None,
        );
    } else {
        for r in no_effect {
            ctx.add(
                NW_NO_EFFECT,
                Severity::Info,
                Element::Reaction(r.to_string()),
                "The reaction changes no variable species: its species are all boundary or \
                 constant, or appear equally on both sides. It only carries a flux.",
                None,
            );
        }
    }
}

fn net_stoichiometry(r: &Reaction) -> BTreeMap<String, f64> {
    let mut net: BTreeMap<String, f64> = BTreeMap::new();
    for (refs, sign) in [(&r.reactants, -1.0), (&r.products, 1.0)] {
        for sr in refs {
            *net.entry(sr.species.clone()).or_default() += sign * sr.stoichiometry.unwrap_or(1.0);
        }
    }
    net.retain(|_, n| n.abs() > 1e-12);
    net
}

/// Every id any expression in the model refers to (rate laws, rules,
/// events, initial assignments, constraints), through function calls.
fn used_symbols(view: &View) -> HashSet<String> {
    let m = view.model;
    let mut used = HashSet::new();
    let mut add = |math: &sbml_io::Math| {
        for s in free_symbols(&view.inline_functions(math)) {
            used.insert(s);
        }
    };
    for r in &m.reactions {
        if let Some(x) = r.kinetic_law.as_ref().and_then(|k| k.math.as_ref()) {
            add(x);
        }
    }
    for rule in &m.rules {
        if let Some(x) = rule.math() {
            add(x);
        }
    }
    for ia in &m.initial_assignments {
        if let Some(x) = &ia.math {
            add(x);
        }
    }
    for e in &m.events {
        for x in e
            .trigger
            .iter()
            .filter_map(|t| t.math.as_ref())
            .chain(e.delay.iter().filter_map(|d| d.math.as_ref()))
            .chain(e.assignments.iter().filter_map(|a| a.math.as_ref()))
        {
            add(x);
        }
    }
    for c in &m.constraints {
        if let Some(x) = &c.math {
            add(x);
        }
    }
    used
}

/// Per-species findings, or one model-level finding listing them all when
/// there are many.
fn report_species(
    ctx: &mut Ctx,
    code: &'static str,
    species: Vec<String>,
    one: impl Fn(&str) -> String,
    many: impl Fn(&str) -> String,
) {
    if species.len() >= GROUP_AT {
        ctx.add(
            code,
            Severity::Info,
            Element::Model,
            many(&species.join(", ")),
            None,
        );
    } else {
        for s in species {
            let msg = one(&s);
            ctx.add(code, Severity::Info, Element::Species(s), msg, None);
        }
    }
}

// ------------------------------------------------------- disconnected

fn find(parent: &mut [usize], i: usize) -> usize {
    let mut r = i;
    while parent[r] != r {
        r = parent[r];
    }
    let mut x = i;
    while parent[x] != r {
        let next = parent[x];
        parent[x] = r;
        x = next;
    }
    r
}

/// NW004: parts of the network that never affect each other.
fn disconnected(ctx: &mut Ctx) {
    let view = ctx.view;
    let m = view.model;
    if m.reactions.len() < 2 {
        return;
    }
    // Nodes: species, then reactions.
    let index: HashMap<&str, usize> = m
        .species
        .iter()
        .enumerate()
        .map(|(i, s)| (s.id.as_str(), i))
        .collect();
    let ns = m.species.len();
    let mut parent: Vec<usize> = (0..ns + m.reactions.len()).collect();
    let union = |a: usize, b: usize, parent: &mut Vec<usize>| {
        let (ra, rb) = (find(parent, a), find(parent, b));
        if ra != rb {
            parent[ra] = rb;
        }
    };
    let species_in = |math: &sbml_io::Math| -> Vec<usize> {
        free_symbols(&view.expand(math, &HashSet::new()))
            .iter()
            .filter_map(|s| index.get(s.as_str()).copied())
            .collect()
    };
    for (ri, r) in m.reactions.iter().enumerate() {
        let node = ns + ri;
        for sr in r.reactants.iter().chain(&r.products) {
            if let Some(&i) = index.get(sr.species.as_str()) {
                union(node, i, &mut parent);
            }
        }
        for md in &r.modifiers {
            if let Some(&i) = index.get(md.species.as_str()) {
                union(node, i, &mut parent);
            }
        }
        if let Some(x) = r.kinetic_law.as_ref().and_then(|k| k.math.as_ref()) {
            for i in species_in(x) {
                union(node, i, &mut parent);
            }
        }
    }
    // Rules and events tie together the species they involve.
    let tie = |ids: Vec<usize>, parent: &mut Vec<usize>| {
        for w in ids.windows(2) {
            union(w[0], w[1], parent);
        }
    };
    for rule in &m.rules {
        let mut ids: Vec<usize> = rule.math().map(species_in).unwrap_or_default();
        if let Some(v) = rule.variable().and_then(|v| index.get(v).copied()) {
            ids.push(v);
        }
        tie(ids, &mut parent);
    }
    for e in &m.events {
        let mut ids = Vec::new();
        for x in e
            .trigger
            .iter()
            .filter_map(|t| t.math.as_ref())
            .chain(e.assignments.iter().filter_map(|a| a.math.as_ref()))
        {
            ids.extend(species_in(x));
        }
        ids.extend(
            e.assignments
                .iter()
                .filter_map(|a| index.get(a.variable.as_str()).copied()),
        );
        tie(ids, &mut parent);
    }
    // Parts with at least one reaction.
    let mut parts: BTreeMap<usize, (Vec<&str>, Vec<&str>)> = BTreeMap::new();
    for (ri, r) in m.reactions.iter().enumerate() {
        let root = find(&mut parent, ns + ri);
        parts.entry(root).or_default().1.push(r.id.as_str());
    }
    if parts.len() < 2 {
        return;
    }
    for (i, s) in m.species.iter().enumerate() {
        let root = find(&mut parent, i);
        if let Some(p) = parts.get_mut(&root) {
            p.0.push(s.id.as_str());
        }
    }
    let mut parts: Vec<(Vec<&str>, Vec<&str>)> = parts.into_values().collect();
    parts.sort_by_key(|p| std::cmp::Reverse(p.1.len()));
    let shown: Vec<String> = parts
        .iter()
        .take(6)
        .enumerate()
        .map(|(i, (sp, rs))| {
            format!(
                "({}) {} reaction{} ({}) with species {}",
                i + 1,
                rs.len(),
                if rs.len() == 1 { "" } else { "s" },
                abbreviate(rs, 4),
                abbreviate(sp, 6)
            )
        })
        .collect();
    let more = parts.len().saturating_sub(6);
    ctx.add(
        NW_DISCONNECTED,
        Severity::Info,
        Element::Model,
        format!(
            "The reaction network falls into {} separate parts that never affect each other: {}{}. \
             This can be intended (independent modules), but is sometimes a missing link or a \
             misspelt species.",
            parts.len(),
            shown.join("; "),
            if more > 0 {
                format!("; and {more} more")
            } else {
                String::new()
            }
        ),
        None,
    );
}

fn abbreviate(items: &[&str], n: usize) -> String {
    if items.len() <= n {
        items.join(", ")
    } else {
        format!("{}, ... ({} in all)", items[..n].join(", "), items.len())
    }
}

// ------------------------------------------------------ conservation

/// The species only reactions change, with their initial amounts.
fn reaction_species(
    view: &View,
    held: &dyn Fn(&str) -> bool,
    set_otherwise: &dyn Fn(&str) -> bool,
) -> Vec<(String, f64, f64)> {
    let m = view.model;
    m.species
        .iter()
        .filter(|s| {
            view.in_reactions.contains(s.id.as_str()) && !held(&s.id) && !set_otherwise(&s.id)
        })
        .map(|s| {
            let x0 = view.range(&s.id).map(|r| r.nominal).unwrap_or(f64::NAN);
            let size = m
                .compartments
                .iter()
                .find(|c| c.id == s.compartment)
                .and_then(|c| view.range(&c.id).map(|r| r.nominal))
                .unwrap_or(1.0);
            // Conservation is of amounts; a concentration times its size.
            let scale = if s.has_only_substance_units == Some(true) {
                1.0
            } else {
                size
            };
            (s.id.clone(), x0, scale)
        })
        .collect()
}

/// NW005: conserved totals.
fn conservation(
    ctx: &mut Ctx,
    nets: &[BTreeMap<String, f64>],
    held: &dyn Fn(&str) -> bool,
    set_otherwise: &dyn Fn(&str) -> bool,
) {
    let view = ctx.view;
    let species = reaction_species(view, held, set_otherwise);
    if species.len() < 2 {
        return;
    }
    // A: one row per reaction, one column per species; c with A c = 0 is a
    // conservation law c^T N = 0 (in amounts).
    let a: Vec<Vec<f64>> = nets
        .iter()
        .map(|n| {
            species
                .iter()
                .map(|(s, _, _)| n.get(s).copied().unwrap_or(0.0))
                .collect()
        })
        .collect();
    let laws = null_space(a, species.len());
    // Tidying compares every pair of laws: fine for a few dozen, too slow
    // for a genome-scale model's hundreds.
    let laws = if laws.len() <= 60 { tidy(laws) } else { laws };
    if laws.is_empty() {
        return;
    }
    let single_size = species
        .iter()
        .map(|s| s.2)
        .all(|v| (v - species[0].2).abs() < 1e-12);
    let mut lines = Vec::new();
    for c in &laws {
        let c = integerize(c);
        // Mostly positive reads best.
        let negatives = c.iter().filter(|x| **x < 0.0).count();
        let positives = c.iter().filter(|x| **x > 0.0).count();
        let c: Vec<f64> = if negatives > positives {
            c.iter().map(|x| -x).collect()
        } else {
            c
        };
        let mut text = String::new();
        let mut total = 0.0;
        let mut first = true;
        for (ci, (s, x0, scale)) in c.iter().zip(&species) {
            if ci.abs() < 1e-9 {
                continue;
            }
            let mag = ci.abs();
            let coef = if (mag - 1.0).abs() < 1e-9 {
                String::new()
            } else {
                format!("{}*", fmt_num(mag))
            };
            let sign = match (first, *ci < 0.0) {
                (true, true) => "-",
                (true, false) => "",
                (false, true) => " - ",
                (false, false) => " + ",
            };
            text.push_str(&format!("{sign}{coef}{s}"));
            first = false;
            // Concentrations in one compartment: the total in those terms.
            let x = if single_size { *x0 } else { x0 * scale };
            total += ci * x;
        }
        let total = if total.is_finite() {
            format!(" = {}", fmt_num(total))
        } else {
            String::new()
        };
        lines.push(format!("{text}{total}"));
    }
    let n = lines.len();
    let shown = lines
        .iter()
        .take(10)
        .cloned()
        .collect::<Vec<_>>()
        .join("\n");
    let more = n.saturating_sub(10);
    ctx.add(
        NW_CONSERVED,
        Severity::Info,
        Element::Model,
        format!(
            "{n} conserved total{} (conservation law{}{}): the reactions move material between \
             these species but never change the total{}.\n{shown}{}",
            if n == 1 { "" } else { "s" },
            if n == 1 { "" } else { "s" },
            if single_size { "" } else { ", in amounts" },
            if single_size {
                ""
            } else {
                "; each species' amount is its concentration times its compartment's size"
            },
            if more > 0 {
                format!("\nand {more} more")
            } else {
                String::new()
            }
        ),
        None,
    );
}

/// A basis of the null space of `a` (rows by `cols`), from its reduced row
/// echelon form.
fn null_space(mut a: Vec<Vec<f64>>, cols: usize) -> Vec<Vec<f64>> {
    let mut pivots: Vec<(usize, usize)> = Vec::new();
    let mut row = 0;
    for col in 0..cols {
        if row >= a.len() {
            break;
        }
        let (best, val) = (row..a.len())
            .map(|r| (r, a[r][col].abs()))
            .fold((row, 0.0), |b, x| if x.1 > b.1 { x } else { b });
        if val < 1e-9 {
            continue;
        }
        a.swap(row, best);
        let p = a[row][col];
        for v in a[row].iter_mut() {
            *v /= p;
        }
        for r in 0..a.len() {
            if r != row && a[r][col].abs() > 1e-12 {
                let f = a[r][col];
                let pivot_row = a[row].clone();
                for (v, pv) in a[r].iter_mut().zip(&pivot_row) {
                    *v -= f * pv;
                }
            }
        }
        pivots.push((row, col));
        row += 1;
    }
    let pivot_cols: Vec<usize> = pivots.iter().map(|&(_, c)| c).collect();
    (0..cols)
        .filter(|c| !pivot_cols.contains(c))
        .map(|free| {
            let mut c = vec![0.0; cols];
            c[free] = 1.0;
            for &(r, pc) in &pivots {
                c[pc] = -a[r][free];
            }
            c
        })
        .collect()
}

/// How readable a conservation law is: fewer negative coefficients (after
/// choosing the better overall sign), then fewer species.
fn badness(v: &[f64]) -> (usize, usize) {
    let neg = v.iter().filter(|x| **x < -1e-9).count();
    let pos = v.iter().filter(|x| **x > 1e-9).count();
    (neg.min(pos), neg + pos)
}

/// Recombines a null-space basis into laws that read naturally: adding or
/// subtracting small multiples of one law to another, keeping a change
/// whenever it removes negative coefficients (or species). `-E + S + P`
/// plus `E + ES` becomes `S + ES + P`.
fn tidy(mut basis: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    for _ in 0..4 {
        let mut improved = false;
        for i in 0..basis.len() {
            for j in 0..basis.len() {
                if i == j {
                    continue;
                }
                for t in [1.0, -1.0, 2.0, -2.0] {
                    let cand: Vec<f64> = basis[i]
                        .iter()
                        .zip(&basis[j])
                        .map(|(a, b)| a + t * b)
                        .collect();
                    if cand.iter().all(|x| x.abs() < 1e-9) {
                        continue;
                    }
                    if badness(&cand) < badness(&basis[i]) {
                        basis[i] = cand;
                        improved = true;
                    }
                }
            }
        }
        if !improved {
            break;
        }
    }
    basis
}

/// Scales a vector to small whole numbers if it is close to one.
fn integerize(c: &[f64]) -> Vec<f64> {
    for d in 1..=12u32 {
        let d = d as f64;
        if c.iter().all(|x| ((x * d) - (x * d).round()).abs() < 1e-7) {
            return c.iter().map(|x| (x * d).round()).collect();
        }
    }
    c.to_vec()
}

// ---------------------------------------------------- never present

/// NW006: species that start at zero and can never be produced, because
/// every reaction that would make one of them needs another as a reactant.
fn never_present(
    ctx: &mut Ctx,
    nets: &[BTreeMap<String, f64>],
    held: &dyn Fn(&str) -> bool,
    set_otherwise: &dyn Fn(&str) -> bool,
) {
    let view = ctx.view;
    let m = view.model;
    let zero = |id: &str| view.range(id).is_some_and(|r| r.nominal == 0.0);
    // Species stuck at zero for good: candidates that start at zero and are
    // changed only by reactions, plus fixed species held at zero.
    let mut stuck: BTreeSet<String> = m
        .species
        .iter()
        .filter(|s| {
            view.in_reactions.contains(s.id.as_str()) && !set_otherwise(&s.id) && zero(&s.id)
        })
        .map(|s| s.id.clone())
        .collect();
    let fixed_zero: HashSet<String> = stuck.iter().filter(|s| held(s)).cloned().collect();
    let laws: Vec<Option<sbml_io::Math>> = m
        .reactions
        .iter()
        .map(|r| super::expanded_rate(view, r))
        .collect();
    // Which way a reaction can run with every stuck species at zero, asked
    // of its rate law (over all other non-negative states): forward
    // (making its products), backward (making its reactants), or not at
    // all. A rate law that ignores a reactant still runs without it.
    let can_run = |ri: usize, stuck: &BTreeSet<String>| -> (bool, bool) {
        let Some(law) = &laws[ri] else {
            return (false, false);
        };
        let mut scope = super::scope::Scope::new(view);
        scope.locals = super::locals(&m.reactions[ri]);
        for s in stuck {
            scope.pins.insert(s.clone(), 0.0);
        }
        let (iv, _) = scope.interval(law);
        (iv.hi > 0.0, iv.lo < 0.0)
    };
    loop {
        let mut changed = false;
        for (ri, net) in nets.iter().enumerate() {
            let (forward, backward) = can_run(ri, &stuck);
            for (s, n) in net {
                let made = (*n > 0.0 && forward) || (*n < 0.0 && backward);
                if made && !fixed_zero.contains(s) && stuck.remove(s) {
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }
    // Report the variable ones; species fixed at zero are by design.
    let final_stuck = stuck.clone();
    let stuck: Vec<String> = stuck.into_iter().filter(|s| !held(s)).collect();
    if stuck.is_empty() {
        return;
    }
    // Reactions involving them that can never run at all.
    let dead: Vec<&str> = m
        .reactions
        .iter()
        .enumerate()
        .filter(|(ri, r)| {
            r.reactants
                .iter()
                .chain(&r.products)
                .any(|sr| stuck.contains(&sr.species))
                && can_run(*ri, &final_stuck) == (false, false)
        })
        .map(|(_, r)| r.id.as_str())
        .collect();
    let names: Vec<&str> = stuck.iter().map(String::as_str).collect();
    ctx.add(
        NW_NEVER_PRESENT,
        Severity::Warning,
        Element::Model,
        format!(
            "{} start{} at zero and can never be made: every reaction that could produce {} has \
             a rate of zero while {} {} zero, so {} stay{} at zero for the whole simulation{}. \
             This may be intended (a knockout, a pathway switched off); if not, an initial \
             value is probably missing, or a reaction that should start {} off.",
            abbreviate(&names, 12),
            if names.len() == 1 { "s" } else { "" },
            if names.len() == 1 {
                "it"
            } else {
                "one of them"
            },
            if names.len() == 1 { "it" } else { "they" },
            if names.len() == 1 { "is" } else { "are" },
            if names.len() == 1 { "it" } else { "they" },
            if names.len() == 1 { "s" } else { "" },
            if dead.is_empty() {
                String::new()
            } else {
                format!(
                    ", and {} ({}) never run{}",
                    if dead.len() == 1 {
                        "the reaction involving it"
                    } else {
                        "the reactions involving them"
                    },
                    abbreviate(&dead, 10),
                    if dead.len() == 1 { "s" } else { "" }
                )
            },
            if names.len() == 1 { "it" } else { "them" }
        ),
        None,
    );
}
