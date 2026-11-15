//! Analysing one expression: its free variables, its range over all
//! reachable states, and a search for concrete bad points.

use std::collections::HashMap;

use sbml_io::{Math, Symbol};

use crate::eval::{Eval, Hazard, TIME};
use crate::sample::{self, Var};
use crate::value::{Interval, Status};
use crate::view::{Range, View};

/// How many concrete points to try per question.
pub const BUDGET: usize = 300;

/// Where an expression is evaluated: the model, plus a rate law's local
/// parameters, plus values pinned for the question being asked.
pub struct Scope<'v, 'm> {
    pub view: &'v View<'m>,
    pub locals: HashMap<String, f64>,
    pub pins: HashMap<String, f64>,
}

impl<'v, 'm> Scope<'v, 'm> {
    pub fn new(view: &'v View<'m>) -> Self {
        Scope {
            view,
            locals: HashMap::new(),
            pins: HashMap::new(),
        }
    }

    pub fn range(&self, id: &str) -> Option<Range> {
        let fixed = |x: f64| Range {
            interval: Interval::point(x),
            nominal: x,
            varies: false,
        };
        if let Some(x) = self.pins.get(id) {
            return Some(fixed(*x));
        }
        if let Some(x) = self.locals.get(id) {
            return Some(if x.is_finite() {
                fixed(*x)
            } else {
                Range {
                    interval: Interval::ALL,
                    nominal: 1.0,
                    varies: true,
                }
            });
        }
        self.view.range(id)
    }

    /// The symbols `m` depends on that vary, in a stable order.
    pub fn vars(&self, m: &Math) -> Vec<Var> {
        let mut ids = free_symbols(m);
        ids.sort();
        ids.dedup();
        ids.into_iter()
            .filter_map(|id| {
                let r = self.range(&id)?;
                r.varies.then_some(Var {
                    id,
                    interval: r.interval,
                    nominal: r.nominal,
                })
            })
            .collect()
    }

    /// Values worth trying: the fixed values and numbers in `m`.
    pub fn extras(&self, m: &Math) -> Vec<f64> {
        let mut out = Vec::new();
        collect_numbers(m, &mut out);
        for id in free_symbols(m) {
            if let Some(r) = self.range(&id)
                && !r.varies
                && r.nominal.is_finite()
            {
                out.push(r.nominal);
            }
        }
        out.retain(|x| x.is_finite() && *x >= 0.0);
        out.sort_by(f64::total_cmp);
        out.dedup();
        out.truncate(24);
        out
    }

    pub fn interval(&self, m: &Math) -> (Interval, Vec<Hazard>) {
        let env = |id: &str| self.range(id).map(|r| r.interval);
        let mut ev = Eval::new(&env);
        let v = ev.eval(m);
        (v, ev.hazards)
    }

    /// Whether a point could be a state of the simulation: the species it
    /// sets to zero must not include a whole conserved moiety.
    pub fn reachable(&self, vars: &[Var], p: &[f64]) -> bool {
        let zeros = vars
            .iter()
            .zip(p)
            .filter(|(_, x)| **x == 0.0)
            .map(|(v, _)| v.id.as_str())
            .chain(
                self.pins
                    .iter()
                    .filter(|(_, x)| **x == 0.0)
                    .map(|(k, _)| k.as_str()),
            );
        !self.view.conservation.cannot_all_be_zero(zeros)
    }

    pub fn at(&self, m: &Math, vars: &[Var], p: &[f64]) -> (f64, Vec<Hazard>) {
        let env = |id: &str| -> Option<f64> {
            if let Some(i) = vars.iter().position(|v| v.id == id) {
                return Some(p[i]);
            }
            self.range(id).map(|r| r.nominal)
        };
        let mut ev = Eval::new(&env);
        let v = ev.eval(m);
        (v, ev.hazards)
    }
}

/// Every identifier `m` depends on, with `time` as [`TIME`].
pub fn free_symbols(m: &Math) -> Vec<String> {
    let mut out = Vec::new();
    m.visit_idents(&mut |id| out.push(id.to_string()));
    visit_symbols(m, &mut |s| {
        if s.url == Symbol::TIME {
            out.push(TIME.to_string());
        }
    });
    out
}

fn visit_symbols(m: &Math, f: &mut impl FnMut(&Symbol)) {
    match m {
        Math::Symbol(s) => f(s),
        Math::Apply(_, a) | Math::Call(_, a) => a.iter().for_each(|x| visit_symbols(x, f)),
        Math::CallSymbol(s, a) => {
            f(s);
            a.iter().for_each(|x| visit_symbols(x, f));
        }
        Math::Root { degree: q, arg } | Math::Log { base: q, arg } => {
            if let Some(q) = q {
                visit_symbols(q, f);
            }
            visit_symbols(arg, f);
        }
        Math::Piecewise { pieces, otherwise } => {
            for (v, c) in pieces {
                visit_symbols(v, f);
                visit_symbols(c, f);
            }
            if let Some(o) = otherwise {
                visit_symbols(o, f);
            }
        }
        Math::Lambda { body, .. } => visit_symbols(body, f),
        Math::Semantics { math, .. } => visit_symbols(math, f),
        _ => {}
    }
}

fn collect_numbers(m: &Math, out: &mut Vec<f64>) {
    match m {
        Math::Number { value, .. } => out.push(value.value()),
        Math::Apply(_, a) | Math::Call(_, a) | Math::CallSymbol(_, a) => {
            a.iter().for_each(|x| collect_numbers(x, out))
        }
        Math::Root { degree: q, arg } | Math::Log { base: q, arg } => {
            if let Some(q) = q {
                collect_numbers(q, out);
            }
            collect_numbers(arg, out);
        }
        Math::Piecewise { pieces, otherwise } => {
            for (v, c) in pieces {
                collect_numbers(v, out);
                collect_numbers(c, out);
            }
            if let Some(o) = otherwise {
                collect_numbers(o, out);
            }
        }
        Math::Lambda { body, .. } => collect_numbers(body, out),
        Math::Semantics { math, .. } => collect_numbers(math, out),
        _ => {}
    }
}

/// A hazard found at a concrete point, with the point as text.
pub struct Witnessed {
    pub hazard: Hazard,
    pub witness: String,
    /// Found at the model's initial values: a simulation fails at once.
    pub at_initial: bool,
}

/// Searches for points where `m` certainly fails. One per distinct
/// hazard, with the first point that showed it.
pub fn search_hazards(scope: &Scope, m: &Math) -> Vec<Witnessed> {
    let vars = scope.vars(m);
    let extras = scope.extras(m);
    let mut found: Vec<Witnessed> = Vec::new();
    let initial: Vec<f64> = vars.iter().map(|v| v.nominal).collect();
    let (_, at_initial) = scope.at(m, &vars, &initial);
    for h in at_initial {
        if h.status == Status::Certain
            && !found
                .iter()
                .any(|w| w.hazard.kind == h.kind && w.hazard.operand == h.operand)
        {
            // The values that matter: those of the offending operand.
            let mut ids = free_symbols(&h.operand_math);
            ids.sort();
            ids.dedup();
            let shown: Vec<String> = ids
                .iter()
                .filter_map(|id| {
                    let r = scope.range(id)?;
                    Some(format!(
                        "{} = {}",
                        sample::display_id(id),
                        sample::fmt_num(r.nominal)
                    ))
                })
                .collect();
            let witness = if shown.is_empty() {
                "at the model's values".to_string()
            } else {
                format!("initial values: {}", shown.join(", "))
            };
            found.push(Witnessed {
                hazard: h,
                witness,
                at_initial: true,
            });
        }
    }
    for p in sample::points(&vars, &extras, BUDGET) {
        if !scope.reachable(&vars, &p) {
            continue;
        }
        let (_, hz) = scope.at(m, &vars, &p);
        for h in hz {
            if h.status != Status::Certain {
                continue;
            }
            if found
                .iter()
                .any(|w| w.hazard.kind == h.kind && w.hazard.operand == h.operand)
            {
                continue;
            }
            let witness = if vars.is_empty() {
                "at the model's values".to_string()
            } else {
                sample::describe(&vars, &p)
            };
            found.push(Witnessed {
                hazard: h,
                witness,
                at_initial: false,
            });
        }
    }
    found
}
