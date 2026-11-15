//! Thermodynamic consistency: TD001 Haldane, TD002 Wegscheider.
//!
//! A reversible rate law is zero at equilibrium, and thermodynamics says
//! that happens at one value of the mass-action ratio
//! `Gamma = prod(x_i ^ n_i)` (net stoichiometry `n_i`, products positive):
//! the equilibrium constant. Catalysts and effectors change how fast the
//! reaction gets there, never where "there" is.
//!
//! TD001: for each rate law that changes sign, the ratio at which it is zero
//! is found numerically (scaling the products until the rate vanishes) at
//! several random states. If the ratio moves, something other than the
//! ratio shifts the equilibrium; perturbing each variable in turn says
//! which. This needs no knowledge of the rate law's form.
//!
//! TD002: where reversible reactions with an equilibrium constant form a
//! closed cycle (a combination returning every species, boundary species
//! included, to where it started), their constants raised to the cycle's
//! coefficients must multiply to 1; otherwise the cycle would turn forever
//! at equilibrium. Cycles come from the null space of the stoichiometry.

use std::collections::BTreeMap;

use sbml_io::Math;

use super::scope::{Scope, free_symbols};
use super::{Ctx, expanded_rate, locals};
use crate::codes::{TD_HALDANE, TD_WEGSCHEIDER};
use crate::eval::TIME;
use crate::finding::{Element, Severity};
use crate::fmt::infix;
use crate::sample::{Var, fmt_num};
use crate::value::Interval;

/// Random states per reaction.
const POINTS: usize = 6;
/// Search range for ln of the products' scale factor.
const LN_RANGE: f64 = 40.0;

pub fn check(ctx: &mut Ctx) {
    let view = ctx.view;
    // Reactions with a consistent equilibrium constant: (index, ln Keq).
    let mut keq: Vec<(usize, f64)> = Vec::new();
    let mut rng = Rng(0x2545_F491_4F6C_DD1D);

    for (ri, r) in view.model.reactions.iter().enumerate() {
        let Some(law) = expanded_rate(view, r) else { continue };
        let deps = free_symbols(&law);
        if deps.iter().any(|d| d == TIME) {
            continue; // an equilibrium that moves with time is not a constant
        }
        let net = net_stoichiometry(r);
        // Only chemical reactions with both sides, whole-number
        // stoichiometries (fractions and large numbers usually encode
        // volume conversions), and a rate law that depends on both sides.
        // One-sided "reactions" (currents, sources) have no equilibrium;
        // a rate law ignoring one side is RL005/RL006's business.
        let elementary = net.values().all(|n| n.fract() == 0.0 && n.abs() <= 4.0);
        let uses = |positive: bool| {
            net.iter()
                .any(|(s, n)| (*n > 0.0) == positive && deps.iter().any(|d| d == s))
        };
        if !elementary || !uses(true) || !uses(false) {
            continue;
        }
        let mut scope = Scope::new(view);
        scope.locals = locals(r);

        // The variables to set: the participants and every other species
        // the rate law uses, even ones the model holds fixed, since whether
        // an enzyme shifts the equilibrium is a property of the rate law.
        // Parameters stay at their values: one changed by an event changes
        // the conditions, not the rate law.
        let mut vars: Vec<Var> = Vec::new();
        for s in net.keys() {
            let nominal = scope.range(s).map(|x| x.nominal).unwrap_or(1.0);
            vars.push(Var {
                id: s.clone(),
                interval: Interval::NON_NEGATIVE,
                nominal,
            });
        }
        for d in &deps {
            if view.kinds.get(d.as_str()) == Some(&crate::view::Kind::Species)
                && !vars.iter().any(|x| &x.id == d)
            {
                let nominal = scope.range(d).map(|x| x.nominal).unwrap_or(1.0);
                vars.push(Var {
                    id: d.clone(),
                    interval: Interval::NON_NEGATIVE,
                    nominal,
                });
            }
        }
        // Scale the products; with none, the reactants.
        let products: Vec<usize> = (0..vars.len())
            .filter(|&i| net.get(&vars[i].id).is_some_and(|n| *n > 0.0))
            .collect();
        let scaled: Vec<usize> = if products.is_empty() {
            (0..vars.len())
                .filter(|&i| net.get(&vars[i].id).is_some_and(|n| *n < 0.0))
                .collect()
        } else {
            products
        };

        let mut values: Vec<(Vec<f64>, f64)> = Vec::new();
        for _ in 0..POINTS {
            let p: Vec<f64> = vars
                .iter()
                .map(|v| {
                    let scale = if v.nominal.is_finite() && v.nominal > 0.0 { v.nominal } else { 1.0 };
                    scale * 10f64.powf(rng.unit() * 2.0 - 1.0)
                })
                .collect();
            if let Some(g) = ln_ratio_at_equilibrium(&scope, &law, &vars, &p, &net, &scaled) {
                values.push((p, g));
            }
        }
        if values.len() < 2 {
            continue; // never zero: not reversible in practice
        }
        let lo = values.iter().map(|x| x.1).fold(f64::INFINITY, f64::min);
        let hi = values.iter().map(|x| x.1).fold(f64::NEG_INFINITY, f64::max);
        if hi - lo <= tolerance(lo, hi) {
            keq.push((ri, (lo + hi) / 2.0));
            continue;
        }

        // The equilibrium moves: find what moves it.
        let (p0, g0) = &values[0];
        let mut culprits: Vec<String> = Vec::new();
        for (i, v) in vars.iter().enumerate() {
            let mut p = p0.clone();
            p[i] *= 3.0;
            if let Some(g) = ln_ratio_at_equilibrium(&scope, &law, &vars, &p, &net, &scaled)
                && (g - g0).abs() > tolerance(g, *g0)
            {
                culprits.push(v.id.clone());
            }
        }
        let ratio = ratio_text(&net);
        let moved_by = if culprits.is_empty() {
            "the state of the model".to_string()
        } else {
            list(&culprits)
        };
        let participants: Vec<&String> = culprits.iter().filter(|c| net.contains_key(*c)).collect();
        let why = if !participants.is_empty() && participants.len() == culprits.len() {
            " The reaction's own species shift it, so the powers in the rate law probably \
             do not match the stoichiometry."
        } else {
            " A modifier (an enzyme, inhibitor or activator) can change how fast a \
             reaction reaches equilibrium, but not where the equilibrium is."
        };
        ctx.add(
            TD_HALDANE,
            Severity::Warning,
            Element::Reaction(r.id.clone()),
            format!(
                "The rate law is reversible, but the mass-action ratio {ratio} at which it is \
                 zero (its equilibrium) is not constant: it ranges from {} to {}, and moves with \
                 {moved_by}. A thermodynamically consistent rate law (one satisfying the Haldane \
                 relationship) has one equilibrium constant.{why}\nRate law: {}",
                fmt_num(lo.exp()),
                fmt_num(hi.exp()),
                infix(&law)
            ),
            None,
        );
    }

    wegscheider(ctx, &keq);
}

/// Relative tolerance on ln values, plus a small absolute floor.
fn tolerance(a: f64, b: f64) -> f64 {
    1e-6 + 1e-6 * a.abs().max(b.abs())
}

fn net_stoichiometry(r: &sbml_io::Reaction) -> BTreeMap<String, f64> {
    let mut net: BTreeMap<String, f64> = BTreeMap::new();
    for (refs, sign) in [(&r.reactants, -1.0), (&r.products, 1.0)] {
        for sr in refs {
            *net.entry(sr.species.clone()).or_default() += sign * sr.stoichiometry.unwrap_or(1.0);
        }
    }
    net.retain(|_, n| n.abs() > 1e-12);
    net
}

/// ln of the mass-action ratio at which the rate is zero, holding all but
/// the scaled species at `p`; `None` if the rate does not change sign.
fn ln_ratio_at_equilibrium(
    scope: &Scope,
    law: &Math,
    vars: &[Var],
    p: &[f64],
    net: &BTreeMap<String, f64>,
    scaled: &[usize],
) -> Option<f64> {
    let rate = |ln_l: f64| -> f64 {
        let f = ln_l.exp();
        let mut q = p.to_vec();
        for &i in scaled {
            q[i] *= f;
        }
        scope.at(law, vars, &q).0
    };
    let (mut a, mut b) = (-LN_RANGE, LN_RANGE);
    let (fa, fb) = (rate(a), rate(b));
    if !fa.is_finite() || !fb.is_finite() || fa == 0.0 || fb == 0.0 || fa.signum() == fb.signum() {
        return None;
    }
    let sa = fa.signum();
    for _ in 0..200 {
        let m = (a + b) / 2.0;
        let fm = rate(m);
        if !fm.is_finite() {
            return None;
        }
        if fm == 0.0 {
            a = m;
            b = m;
            break;
        }
        if fm.signum() == sa {
            a = m;
        } else {
            b = m;
        }
    }
    let ln_l = (a + b) / 2.0;
    // ln Gamma at the root: the scaled species are multiplied by e^ln_l.
    let mut g = 0.0;
    for (i, v) in vars.iter().enumerate() {
        if let Some(n) = net.get(&v.id) {
            let x = if scaled.contains(&i) { p[i] * ln_l.exp() } else { p[i] };
            if x <= 0.0 {
                return None;
            }
            g += n * x.ln();
        }
    }
    Some(g)
}

/// `[P]^2/[S]`, for messages.
fn ratio_text(net: &BTreeMap<String, f64>) -> String {
    let term = |s: &str, n: f64| {
        if (n.abs() - 1.0).abs() < 1e-12 {
            format!("[{s}]")
        } else {
            format!("[{s}]^{}", fmt_num(n.abs()))
        }
    };
    let top: Vec<String> = net.iter().filter(|(_, n)| **n > 0.0).map(|(s, n)| term(s, *n)).collect();
    let bottom: Vec<String> = net.iter().filter(|(_, n)| **n < 0.0).map(|(s, n)| term(s, *n)).collect();
    let top = if top.is_empty() { "1".to_string() } else { top.join("*") };
    match bottom.len() {
        0 => top,
        1 => format!("{top}/{}", bottom[0]),
        _ => format!("{top}/({})", bottom.join("*")),
    }
}

fn list(items: &[String]) -> String {
    match items {
        [] => String::new(),
        [a] => a.clone(),
        [a, b] => format!("{a} and {b}"),
        _ => format!("{} and {}", items[..items.len() - 1].join(", "), items[items.len() - 1]),
    }
}

/// TD002: equilibrium constants around closed cycles of reversible
/// reactions must multiply to 1.
fn wegscheider(ctx: &mut Ctx, keq: &[(usize, f64)]) {
    if keq.len() < 2 {
        return;
    }
    let reactions = &ctx.view.model.reactions;
    // Rows: every species (boundary included); columns: the reactions.
    let nets: Vec<BTreeMap<String, f64>> =
        keq.iter().map(|(i, _)| net_stoichiometry(&reactions[*i])).collect();
    let mut species: Vec<&String> = nets.iter().flat_map(|n| n.keys()).collect();
    species.sort();
    species.dedup();
    let k = keq.len();
    let mut a: Vec<Vec<f64>> = species
        .iter()
        .map(|s| nets.iter().map(|n| n.get(*s).copied().unwrap_or(0.0)).collect())
        .collect();
    // Reduced row echelon form; null-space vectors from the free columns.
    let mut pivots: Vec<(usize, usize)> = Vec::new();
    let mut row = 0;
    for col in 0..k {
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
    // A large network can have many failing cycles: report a few, then say
    // how many more there are.
    const SHOWN: usize = 5;
    let mut failing = 0usize;
    for free in (0..k).filter(|c| !pivot_cols.contains(c)) {
        let mut c = vec![0.0; k];
        c[free] = 1.0;
        for &(r, pc) in &pivots {
            c[pc] = -a[r][free];
        }
        // Scale so the smallest coefficient is 1 in size (often integers).
        let min = c.iter().filter(|x| x.abs() > 1e-9).map(|x| x.abs()).fold(f64::INFINITY, f64::min);
        for x in c.iter_mut() {
            *x /= min;
        }
        let sum: f64 = c.iter().zip(keq).map(|(ci, (_, lk))| ci * lk).sum();
        // Within 1% of 1 is published constants rounded to a few digits,
        // not an inconsistency worth a warning.
        if sum.abs() <= 0.01 {
            continue;
        }
        failing += 1;
        if failing > SHOWN {
            continue;
        }
        let members: Vec<(f64, usize, f64)> = c
            .iter()
            .zip(keq)
            .filter(|(ci, _)| ci.abs() > 1e-9)
            .map(|(ci, (ri, lk))| (*ci, *ri, *lk))
            .collect();
        let cycle = members
            .iter()
            .enumerate()
            .map(|(j, (ci, ri, _))| {
                let id = &reactions[*ri].id;
                let mag = ci.abs();
                let coef = if (mag - 1.0).abs() < 1e-9 { String::new() } else { format!("{} ", fmt_num(mag)) };
                let sign = match (j, *ci < 0.0) {
                    (0, true) => "-".to_string(),
                    (0, false) => String::new(),
                    (_, true) => " - ".to_string(),
                    (_, false) => " + ".to_string(),
                };
                format!("{sign}{coef}{id}")
            })
            .collect::<String>();
        let constants = members
            .iter()
            .map(|(_, ri, lk)| format!("Keq({}) = {}", reactions[*ri].id, fmt_num(lk.exp())))
            .collect::<Vec<_>>()
            .join(", ");
        let names: Vec<String> = members.iter().map(|(_, ri, _)| reactions[*ri].id.clone()).collect();
        ctx.add(
            TD_WEGSCHEIDER,
            Severity::Warning,
            Element::Model,
            format!(
                "The reversible reactions {} form a closed cycle: running {cycle} (a minus sign \
                 means that reaction runs backwards) returns every species to where it started. \
                 Their equilibrium constants, raised to those coefficients, must then multiply \
                 to 1 (the Wegscheider condition), but they multiply to {}: {constants}. Such a \
                 cycle would keep turning at equilibrium, which is thermodynamically impossible.",
                list(&names),
                fmt_num(sum.exp()),
            ),
            None,
        );
    }
    if failing > SHOWN {
        let more = failing - SHOWN;
        ctx.add(
            TD_WEGSCHEIDER,
            Severity::Warning,
            Element::Model,
            format!(
                "{more} more cycle{} of reversible reactions also fail{} the Wegscheider \
                 condition; only the first {SHOWN} are listed.",
                if more == 1 { "" } else { "s" },
                if more == 1 { "s" } else { "" }
            ),
            None,
        );
    }
}

/// A small deterministic random number generator.
struct Rng(u64);

impl Rng {
    fn unit(&mut self) -> f64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        (x >> 11) as f64 / (1u64 << 53) as f64
    }
}
