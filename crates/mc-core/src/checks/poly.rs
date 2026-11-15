//! Rate laws as polynomials in the species: the shape mass action has.

use std::collections::BTreeMap;

use sbml_io::{Math, Op};

use super::scope::free_symbols;

/// A term of a mass-action expression: its sign and species exponents.
pub type Mono = (f64, BTreeMap<String, f64>);

const MAX_TERMS: usize = 32;

/// `m` as a sum of signed monomials in the species, if it is one. Factors
/// free of species (rate constants, compartments) are dropped, keeping
/// their sign.
pub fn monomials(
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
