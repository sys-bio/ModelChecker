//! Concrete points to try when looking for a counterexample.
//!
//! Interval analysis can say a problem is possible; a point where it really
//! happens turns that into a definite finding, with a witness to show. The
//! points mix the obvious corners (everything zero, everything at its
//! initial value) with values chosen to make differences vanish (a species
//! equal to a parameter appearing in the same expression) and a fixed
//! pseudo-random spread.

use crate::value::Interval;

/// A free variable of an expression and what it may range over.
#[derive(Debug, Clone)]
pub struct Var {
    pub id: String,
    pub interval: Interval,
    pub nominal: f64,
}

/// Points to try: each is one value per variable, in the order of `vars`.
/// `extras` are other values worth trying (the expression's constants).
pub fn points(vars: &[Var], extras: &[f64], budget: usize) -> Vec<Vec<f64>> {
    let n = vars.len();
    let mut out: Vec<Vec<f64>> = Vec::new();
    let clamp = |v: &Var, x: f64| -> Option<f64> {
        (x.is_finite() && v.interval.contains(x)).then_some(x)
    };
    let base_nominal: Vec<f64> = vars
        .iter()
        .map(|v| clamp(v, v.nominal).unwrap_or_else(|| low(v)))
        .collect();
    let base_low: Vec<f64> = vars.iter().map(low).collect();
    let push = |p: Vec<f64>, out: &mut Vec<Vec<f64>>| {
        if out.len() < budget && !out.contains(&p) {
            out.push(p);
        }
    };
    push(base_nominal.clone(), &mut out);
    push(base_low.clone(), &mut out);
    if n == 0 {
        return out;
    }

    let mut candidates: Vec<f64> = vec![0.0, 1.0, 0.5, 2.0, 10.0, 1e-3, 1e3];
    for &x in extras {
        if x.is_finite() && !candidates.contains(&x) {
            candidates.push(x);
        }
    }
    // Everything at one value: makes differences like S - P vanish.
    for &c in &candidates {
        let p: Vec<f64> = vars.iter().map(|v| clamp(v, c).unwrap_or(low(v))).collect();
        push(p, &mut out);
    }
    // One variable moved, the others at their initial values or at zero.
    for (i, v) in vars.iter().enumerate() {
        for &c in &candidates {
            let Some(x) = clamp(v, c) else { continue };
            for base in [&base_nominal, &base_low] {
                let mut p = base.clone();
                p[i] = x;
                push(p, &mut out);
            }
        }
    }
    // Pseudo-random spread, the same every run.
    let mut rng = XorShift(0x9E37_79B9_7F4A_7C15);
    while out.len() < budget {
        let before = out.len();
        let p: Vec<f64> = vars
            .iter()
            .map(|v| {
                let r = rng.next();
                let bounded = v.interval.lo.is_finite() && v.interval.hi.is_finite();
                let x = match r % 4 {
                    0 => 0.0,
                    1 if !extras.is_empty() => extras[(r >> 8) as usize % extras.len()],
                    // Spread evenly through a bounded range...
                    2 if bounded => {
                        v.interval.lo + rng.unit() * (v.interval.hi - v.interval.lo)
                    }
                    // ...and over orders of magnitude otherwise.
                    _ => 10f64.powf(rng.unit() * 8.0 - 4.0),
                };
                clamp(v, x).unwrap_or_else(|| clamp(v, -x).unwrap_or(low(v)))
            })
            .collect();
        push(p, &mut out);
        if out.len() == before {
            // Few distinct points exist; stop rather than spin.
            if rng.0.is_multiple_of(64) {
                break;
            }
        }
    }
    out
}

fn low(v: &Var) -> f64 {
    if v.interval.lo.is_finite() {
        v.interval.lo
    } else if v.interval.contains(0.0) {
        0.0
    } else {
        v.interval.hi
    }
}

struct XorShift(u64);

impl XorShift {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Text for a point, for a witness: `S1 = 0, P = 1`.
pub fn describe(vars: &[Var], p: &[f64]) -> String {
    vars.iter()
        .zip(p)
        .map(|(v, x)| format!("{} = {}", display_id(&v.id), fmt_num(*x)))
        .collect::<Vec<_>>()
        .join(", ")
}

pub fn display_id(id: &str) -> &str {
    match id {
        crate::eval::TIME => "time",
        crate::eval::AVOGADRO => "avogadro",
        _ => id,
    }
}

pub fn fmt_num(x: f64) -> String {
    if x == 0.0 {
        "0".into()
    } else if x.fract() == 0.0 && x.abs() < 1e9 {
        format!("{}", x as i64)
    } else if (1e-3..1e6).contains(&x.abs()) {
        let s = format!("{x:.6}");
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    } else {
        format!("{x:.4e}")
    }
}
