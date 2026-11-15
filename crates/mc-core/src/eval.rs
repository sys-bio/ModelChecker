//! Evaluating an expression in a [`Value`] domain, recording every hazard
//! met on the way (division by zero, log of a non-positive number, ...).
//!
//! Evaluation is lazy the way a simulator's is: only the piecewise branch
//! whose condition holds is evaluated. Over intervals a condition may be
//! undecided; then every branch that might run is evaluated, and hazards in
//! them are marked `guarded`, since the condition may be what protects them
//! (`piecewise(0, S == 0, k/S)`).

use sbml_io::{Constant, Math, Op, Symbol};

use crate::fmt::infix;
use crate::value::{HazardKind, Status, Tri, Value};

/// Names the environment is asked for, for SBML's `csymbol`s.
pub const TIME: &str = "\u{0}time";
pub const AVOGADRO: &str = "\u{0}avogadro";

#[derive(Debug, Clone, PartialEq)]
pub struct Hazard {
    pub kind: HazardKind,
    pub status: Status,
    /// The offending operand, as text: the denominator, the log's argument.
    pub operand: String,
    /// The same, as an expression.
    pub operand_math: Math,
    /// Inside a piecewise branch whose condition was undecided.
    pub guarded: bool,
}

pub struct Eval<'e, V> {
    env: &'e dyn Fn(&str) -> Option<V>,
    pub hazards: Vec<Hazard>,
    guard: usize,
}

impl<'e, V: Value> Eval<'e, V> {
    pub fn new(env: &'e dyn Fn(&str) -> Option<V>) -> Self {
        Eval {
            env,
            hazards: Vec::new(),
            guard: 0,
        }
    }

    fn note(&mut self, kind: HazardKind, status: Status, operand: &Math) {
        if status == Status::Safe {
            return;
        }
        let operand_math = operand.clone();
        let operand = infix(operand);
        let guarded = self.guard > 0;
        // One entry per distinct problem, keeping the strongest status.
        if let Some(h) = self
            .hazards
            .iter_mut()
            .find(|h| h.kind == kind && h.operand == operand)
        {
            if status > h.status {
                h.status = status;
            }
            h.guarded &= guarded;
            return;
        }
        self.hazards.push(Hazard {
            kind,
            status,
            operand,
            operand_math,
            guarded,
        });
    }

    pub fn eval(&mut self, m: &Math) -> V {
        match m {
            Math::Number { value, .. } => V::num(value.value()),
            Math::Ident(id) => (self.env)(id).unwrap_or_else(V::unknown),
            Math::Symbol(s) => self.symbol(s),
            Math::Constant(c) => V::num(match c {
                Constant::True => 1.0,
                Constant::False => 0.0,
                Constant::Pi => std::f64::consts::PI,
                Constant::ExponentialE => std::f64::consts::E,
                Constant::NotANumber => return V::unknown(),
                Constant::Infinity => f64::INFINITY,
            }),
            Math::Apply(op, args) => self.apply(*op, args),
            // Function calls are inlined before evaluation; one left here
            // refers to an undefined function.
            Math::Call(..) => V::unknown(),
            Math::CallSymbol(s, args) => {
                if s.url == Symbol::DELAY && !args.is_empty() {
                    // delay(x, d) takes values x has taken.
                    self.eval(&args[0])
                } else {
                    V::unknown()
                }
            }
            Math::Root { degree, arg } => {
                let x = self.eval(arg);
                let d = match degree {
                    Some(d) => self.eval(d),
                    None => V::num(2.0),
                };
                let (inv, _) = V::num(1.0).div(&d);
                let (r, s, kind) = x.pow(&inv);
                let kind = match kind {
                    HazardKind::NegativeBaseFractionalPower => HazardKind::RootOfNegative,
                    k => k,
                };
                self.note(kind, s, arg);
                r
            }
            Math::Log { base, arg } => {
                let x = self.eval(arg);
                let (lx, s) = x.ln();
                self.note(HazardKind::LogOfNonPositive, s, arg);
                let lb = match base {
                    Some(b) => {
                        let bv = self.eval(b);
                        let (lb, s) = bv.ln();
                        self.note(HazardKind::LogOfNonPositive, s, b);
                        lb
                    }
                    None => V::num(std::f64::consts::LN_10),
                };
                lx.div(&lb).0
            }
            Math::Piecewise { pieces, otherwise } => self.piecewise(pieces, otherwise.as_deref()),
            Math::Lambda { body, .. } => self.eval(body),
            Math::Semantics { math, .. } => self.eval(math),
            Math::Raw(_) => V::unknown(),
        }
    }

    fn symbol(&mut self, s: &Symbol) -> V {
        let key = if s.url == Symbol::TIME {
            TIME
        } else if s.url == Symbol::AVOGADRO {
            AVOGADRO
        } else {
            return V::unknown();
        };
        (self.env)(key).unwrap_or_else(V::unknown)
    }

    fn piecewise(&mut self, pieces: &[(Math, Math)], otherwise: Option<&Math>) -> V {
        let mut acc: Option<V> = None;
        let mut undecided = 0;
        let join = |acc: Option<V>, v: V| Some(match acc {
            Some(a) => a.hull(&v),
            None => v,
        });
        for (value, cond) in pieces {
            let t = self.eval(cond).truth();
            match t {
                Tri::False => continue,
                Tri::True => {
                    self.guard += undecided;
                    let v = self.eval(value);
                    self.guard -= undecided;
                    return join(acc, v).unwrap();
                }
                Tri::Unknown => {
                    self.guard += 1;
                    let v = self.eval(value);
                    self.guard -= 1;
                    acc = join(acc, v);
                    undecided += 1;
                }
            }
        }
        match otherwise {
            Some(o) => {
                self.guard += undecided;
                let v = self.eval(o);
                self.guard -= undecided;
                join(acc, v).unwrap()
            }
            None => acc.unwrap_or_else(V::unknown),
        }
    }

    fn apply(&mut self, op: Op, args: &[Math]) -> V {
        match op {
            Op::Plus => {
                let mut acc = V::num(0.0);
                for a in args {
                    acc = acc.add(&self.eval(a));
                }
                acc
            }
            Op::Times => {
                let mut acc = V::num(1.0);
                for a in args {
                    acc = acc.mul(&self.eval(a));
                }
                acc
            }
            Op::Minus => match args {
                [a] => self.eval(a).neg(),
                [a, b] => {
                    let a = self.eval(a);
                    a.sub(&self.eval(b))
                }
                _ => V::unknown(),
            },
            Op::Divide => match args {
                [a, b] => {
                    let x = self.eval(a);
                    let y = self.eval(b);
                    let (r, s) = x.div(&y);
                    self.note(HazardKind::DivideByZero, s, b);
                    r
                }
                _ => V::unknown(),
            },
            Op::Power => match args {
                [a, b] => {
                    let x = self.eval(a);
                    let y = self.eval(b);
                    let (r, s, kind) = x.pow(&y);
                    self.note(kind, s, a);
                    r
                }
                _ => V::unknown(),
            },
            Op::Ln => match args {
                [a] => {
                    let (r, s) = self.eval(a).ln();
                    self.note(HazardKind::LogOfNonPositive, s, a);
                    r
                }
                _ => V::unknown(),
            },
            Op::Min | Op::Max => {
                let mut it = args.iter();
                let Some(first) = it.next() else {
                    return V::unknown();
                };
                let mut acc = self.eval(first);
                for a in it {
                    let v = self.eval(a);
                    acc = if op == Op::Min { acc.min(&v) } else { acc.max(&v) };
                }
                acc
            }
            Op::Eq | Op::Neq | Op::Gt | Op::Lt | Op::Geq | Op::Leq => {
                // Chains like a < b < c mean a < b and b < c.
                let vals: Vec<V> = args.iter().map(|a| self.eval(a)).collect();
                let mut t = Tri::True;
                for w in vals.windows(2) {
                    t = and(t, V::compare(op, &w[0], &w[1]).truth());
                }
                V::from_tri(t)
            }
            Op::And | Op::Or | Op::Xor | Op::Implies => {
                let ts: Vec<Tri> = args.iter().map(|a| self.eval(a).truth()).collect();
                let t = match op {
                    Op::And => ts.iter().fold(Tri::True, |a, b| and(a, *b)),
                    Op::Or => ts.iter().fold(Tri::False, |a, b| or(a, *b)),
                    Op::Xor => ts.iter().fold(Tri::False, |a, b| xor(a, *b)),
                    _ => match ts.as_slice() {
                        [a, b] => or(not(*a), *b),
                        _ => Tri::Unknown,
                    },
                };
                V::from_tri(t)
            }
            Op::Quotient | Op::Rem => match args {
                [a, b] => {
                    let _ = self.eval(a);
                    let y = self.eval(b);
                    let (_, s) = V::num(1.0).div(&y);
                    self.note(HazardKind::DivideByZero, s, b);
                    V::unknown()
                }
                _ => V::unknown(),
            },
            _ => match args {
                [a] => {
                    let x = self.eval(a);
                    let (r, s) = V::func(op, &x);
                    self.note(HazardKind::Domain, s, a);
                    r
                }
                _ => V::unknown(),
            },
        }
    }
}

fn not(a: Tri) -> Tri {
    match a {
        Tri::True => Tri::False,
        Tri::False => Tri::True,
        Tri::Unknown => Tri::Unknown,
    }
}

fn and(a: Tri, b: Tri) -> Tri {
    match (a, b) {
        (Tri::False, _) | (_, Tri::False) => Tri::False,
        (Tri::True, Tri::True) => Tri::True,
        _ => Tri::Unknown,
    }
}

fn or(a: Tri, b: Tri) -> Tri {
    match (a, b) {
        (Tri::True, _) | (_, Tri::True) => Tri::True,
        (Tri::False, Tri::False) => Tri::False,
        _ => Tri::Unknown,
    }
}

fn xor(a: Tri, b: Tri) -> Tri {
    match (a, b) {
        (Tri::Unknown, _) | (_, Tri::Unknown) => Tri::Unknown,
        (x, y) if x == y => Tri::False,
        _ => Tri::True,
    }
}
