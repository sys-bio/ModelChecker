//! The two number domains expressions are evaluated in: plain `f64`, for
//! concrete points, and [`Interval`], for whole ranges of values at once.
//!
//! Each operation that can fail (division, logarithm, power, ...) reports
//! a [`Status`]: whether it is safe, may fail somewhere in the range, or
//! certainly fails.

use sbml_io::Op;

/// What can go wrong when evaluating an operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum HazardKind {
    DivideByZero,
    /// `0^p` with `p < 0`.
    ZeroToNegativePower,
    /// `x^p` with `x < 0` and `p` not an integer.
    NegativeBaseFractionalPower,
    /// `ln`/`log` of a number <= 0.
    LogOfNonPositive,
    /// Square or other even root of a negative number.
    RootOfNegative,
    /// Any other operation giving NaN from ordinary numbers
    /// (`arcsin(2)`, ...).
    Domain,
}

impl HazardKind {
    pub fn describe(self) -> &'static str {
        match self {
            HazardKind::DivideByZero => "division by zero",
            HazardKind::ZeroToNegativePower => "zero raised to a negative power",
            HazardKind::NegativeBaseFractionalPower => {
                "negative number raised to a fractional power"
            }
            HazardKind::LogOfNonPositive => "logarithm of a non-positive number",
            HazardKind::RootOfNegative => "root of a negative number",
            HazardKind::Domain => "function evaluated outside its domain",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Status {
    Safe,
    Maybe,
    Certain,
}

/// Three-valued truth, for conditions over ranges.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tri {
    True,
    False,
    Unknown,
}

/// A number domain.
pub trait Value: Clone + std::fmt::Debug {
    fn num(x: f64) -> Self;
    /// Anything at all: used for what cannot be analysed.
    fn unknown() -> Self;
    fn add(&self, o: &Self) -> Self;
    fn sub(&self, o: &Self) -> Self;
    fn mul(&self, o: &Self) -> Self;
    fn neg(&self) -> Self;
    fn div(&self, o: &Self) -> (Self, Status);
    /// The result, with the status and kind of any hazard.
    fn pow(&self, o: &Self) -> (Self, Status, HazardKind);
    fn ln(&self) -> (Self, Status);
    /// Built-in one-argument functions other than `ln`.
    fn func(op: Op, x: &Self) -> (Self, Status);
    fn min(&self, o: &Self) -> Self;
    fn max(&self, o: &Self) -> Self;
    /// A relational operator, giving a truth value as a number.
    fn compare(op: Op, a: &Self, b: &Self) -> Self;
    fn truth(&self) -> Tri;
    fn from_tri(t: Tri) -> Self;
    /// The smallest value covering both (for branches that may both run).
    fn hull(&self, o: &Self) -> Self;
}

fn is_int(x: f64) -> bool {
    x.is_finite() && x.fract() == 0.0
}

// ----------------------------------------------------------------- f64

impl Value for f64 {
    fn num(x: f64) -> Self {
        x
    }
    fn unknown() -> Self {
        f64::NAN
    }
    fn add(&self, o: &Self) -> Self {
        self + o
    }
    fn sub(&self, o: &Self) -> Self {
        self - o
    }
    fn mul(&self, o: &Self) -> Self {
        self * o
    }
    fn neg(&self) -> Self {
        -self
    }
    fn div(&self, o: &Self) -> (Self, Status) {
        let s = if *o == 0.0 {
            Status::Certain
        } else {
            Status::Safe
        };
        (self / o, s)
    }
    fn pow(&self, o: &Self) -> (Self, Status, HazardKind) {
        let (x, p) = (*self, *o);
        if x == 0.0 && p < 0.0 {
            (f64::INFINITY, Status::Certain, HazardKind::ZeroToNegativePower)
        } else if x < 0.0 && !is_int(p) && p.is_finite() {
            (f64::NAN, Status::Certain, HazardKind::NegativeBaseFractionalPower)
        } else {
            (x.powf(p), Status::Safe, HazardKind::Domain)
        }
    }
    fn ln(&self) -> (Self, Status) {
        let s = if *self <= 0.0 {
            Status::Certain
        } else {
            Status::Safe
        };
        (f64::ln(*self), s)
    }
    fn func(op: Op, x: &Self) -> (Self, Status) {
        let x = *x;
        let r = match op {
            Op::Abs => x.abs(),
            Op::Exp => x.exp(),
            Op::Floor => x.floor(),
            Op::Ceiling => x.ceil(),
            Op::Factorial => {
                if x < 0.0 || !is_int(x) {
                    f64::NAN
                } else {
                    (1..=(x.min(170.0) as u64)).fold(1.0, |a, k| a * k as f64)
                }
            }
            Op::Sin => x.sin(),
            Op::Cos => x.cos(),
            Op::Tan => x.tan(),
            Op::Sec => 1.0 / x.cos(),
            Op::Csc => 1.0 / x.sin(),
            Op::Cot => 1.0 / x.tan(),
            Op::Sinh => x.sinh(),
            Op::Cosh => x.cosh(),
            Op::Tanh => x.tanh(),
            Op::Sech => 1.0 / x.cosh(),
            Op::Csch => 1.0 / x.sinh(),
            Op::Coth => 1.0 / x.tanh(),
            Op::ArcSin => x.asin(),
            Op::ArcCos => x.acos(),
            Op::ArcTan => x.atan(),
            Op::ArcSec => (1.0 / x).acos(),
            Op::ArcCsc => (1.0 / x).asin(),
            Op::ArcCot => (1.0 / x).atan(),
            Op::ArcSinh => x.asinh(),
            Op::ArcCosh => x.acosh(),
            Op::ArcTanh => x.atanh(),
            Op::ArcSech => (1.0 / x).acosh(),
            Op::ArcCsch => (1.0 / x).asinh(),
            Op::ArcCoth => (1.0 / x).atanh(),
            Op::Not => {
                if x != 0.0 {
                    0.0
                } else {
                    1.0
                }
            }
            _ => f64::NAN,
        };
        let s = if r.is_nan() && !x.is_nan() {
            Status::Certain
        } else {
            Status::Safe
        };
        (r, s)
    }
    fn min(&self, o: &Self) -> Self {
        f64::min(*self, *o)
    }
    fn max(&self, o: &Self) -> Self {
        f64::max(*self, *o)
    }
    fn compare(op: Op, a: &Self, b: &Self) -> Self {
        let t = match op {
            Op::Eq => a == b,
            Op::Neq => a != b,
            Op::Gt => a > b,
            Op::Lt => a < b,
            Op::Geq => a >= b,
            Op::Leq => a <= b,
            _ => false,
        };
        if t { 1.0 } else { 0.0 }
    }
    fn truth(&self) -> Tri {
        if self.is_nan() {
            Tri::Unknown
        } else if *self != 0.0 {
            Tri::True
        } else {
            Tri::False
        }
    }
    fn from_tri(t: Tri) -> Self {
        match t {
            Tri::True => 1.0,
            Tri::False => 0.0,
            Tri::Unknown => f64::NAN,
        }
    }
    fn hull(&self, _o: &Self) -> Self {
        *self
    }
}

// ------------------------------------------------------------ Interval

/// A closed interval of reals, possibly unbounded. Never empty.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interval {
    pub lo: f64,
    pub hi: f64,
}

impl Interval {
    pub const ALL: Interval = Interval {
        lo: f64::NEG_INFINITY,
        hi: f64::INFINITY,
    };
    pub const NON_NEGATIVE: Interval = Interval {
        lo: 0.0,
        hi: f64::INFINITY,
    };

    pub fn new(lo: f64, hi: f64) -> Self {
        if lo.is_nan() || hi.is_nan() || lo > hi {
            Interval::ALL
        } else {
            Interval { lo, hi }
        }
    }
    pub fn point(x: f64) -> Self {
        Interval::new(x, x)
    }
    pub fn contains(&self, x: f64) -> bool {
        self.lo <= x && x <= self.hi
    }
    pub fn is_point(&self) -> bool {
        self.lo == self.hi
    }

    fn from_corners(c: [f64; 4]) -> Self {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for x in c {
            if x.is_nan() {
                return Interval::ALL;
            }
            lo = lo.min(x);
            hi = hi.max(x);
        }
        Interval::new(lo, hi)
    }
}

/// Product where `0 * inf = 0`, the right rule for interval endpoints.
fn emul(a: f64, b: f64) -> f64 {
    if a == 0.0 || b == 0.0 { 0.0 } else { a * b }
}

fn monotone_up(x: Interval, f: impl Fn(f64) -> f64) -> Interval {
    Interval::new(f(x.lo), f(x.hi))
}

fn monotone_down(x: Interval, f: impl Fn(f64) -> f64) -> Interval {
    Interval::new(f(x.hi), f(x.lo))
}

impl Value for Interval {
    fn num(x: f64) -> Self {
        Interval::point(x)
    }
    fn unknown() -> Self {
        Interval::ALL
    }
    fn add(&self, o: &Self) -> Self {
        Interval::new(self.lo + o.lo, self.hi + o.hi)
    }
    fn sub(&self, o: &Self) -> Self {
        Interval::new(self.lo - o.hi, self.hi - o.lo)
    }
    fn mul(&self, o: &Self) -> Self {
        Interval::from_corners([
            emul(self.lo, o.lo),
            emul(self.lo, o.hi),
            emul(self.hi, o.lo),
            emul(self.hi, o.hi),
        ])
    }
    fn neg(&self) -> Self {
        Interval::new(-self.hi, -self.lo)
    }
    fn div(&self, o: &Self) -> (Self, Status) {
        if o.lo > 0.0 || o.hi < 0.0 {
            let inv = Interval::new(1.0 / o.hi, 1.0 / o.lo);
            return (self.mul(&inv), Status::Safe);
        }
        if o.lo == 0.0 && o.hi == 0.0 {
            return (Interval::ALL, Status::Certain);
        }
        // The denominator's range includes zero. If it touches zero only
        // at one end and the numerator has a sign, the result is one-sided.
        let r = if o.lo == 0.0 && self.lo > 0.0 {
            Interval::new(self.lo / o.hi, f64::INFINITY)
        } else if o.lo == 0.0 && self.hi < 0.0 {
            Interval::new(f64::NEG_INFINITY, self.hi / o.hi)
        } else if o.hi == 0.0 && self.lo > 0.0 {
            Interval::new(f64::NEG_INFINITY, self.lo / o.lo)
        } else if o.hi == 0.0 && self.hi < 0.0 {
            Interval::new(self.hi / o.lo, f64::INFINITY)
        } else {
            Interval::ALL
        };
        (r, Status::Maybe)
    }
    fn pow(&self, o: &Self) -> (Self, Status, HazardKind) {
        let x = *self;
        if o.is_point() && is_int(o.lo) {
            let n = o.lo;
            if n >= 0.0 {
                let even = (n / 2.0).fract() == 0.0;
                let r = if !even || x.lo >= 0.0 {
                    monotone_up(x, |v| v.powf(n))
                } else if x.hi <= 0.0 {
                    monotone_down(x, |v| v.powf(n))
                } else {
                    Interval::new(0.0, x.lo.abs().max(x.hi).powf(n))
                };
                return (r, Status::Safe, HazardKind::Domain);
            }
            // Negative integer power: 1 / x^|n|.
            let (d, _, _) = x.pow(&Interval::point(-n));
            let (r, s) = Interval::point(1.0).div(&d);
            return (r, s, HazardKind::ZeroToNegativePower);
        }
        // Fractional or uncertain exponent: needs a non-negative base.
        let mut status = Status::Safe;
        let mut kind = HazardKind::Domain;
        let exact_int_possible = o.is_point(); // a point here is non-integer
        if x.hi < 0.0 && exact_int_possible {
            return (
                Interval::ALL,
                Status::Certain,
                HazardKind::NegativeBaseFractionalPower,
            );
        }
        if x.lo < 0.0 {
            return (
                Interval::ALL,
                Status::Maybe,
                HazardKind::NegativeBaseFractionalPower,
            );
        }
        if x.lo == 0.0 && o.lo < 0.0 {
            status = if x.hi == 0.0 && o.hi < 0.0 {
                Status::Certain
            } else {
                Status::Maybe
            };
            kind = HazardKind::ZeroToNegativePower;
        }
        // For x >= 0, x^y = exp(y ln x) is monotone in each argument, so the
        // extremes are at the corners.
        let p = |a: f64, b: f64| {
            if a == 0.0 {
                if b > 0.0 {
                    0.0
                } else if b == 0.0 {
                    1.0
                } else {
                    f64::INFINITY
                }
            } else {
                a.powf(b)
            }
        };
        let r = Interval::from_corners([p(x.lo, o.lo), p(x.lo, o.hi), p(x.hi, o.lo), p(x.hi, o.hi)]);
        // Exponent ranges crossing zero with base crossing 1 can make the
        // corners miss the value 1; include it to stay sound.
        let r = if o.lo < 0.0 && o.hi > 0.0 {
            r.hull(&Interval::point(1.0))
        } else {
            r
        };
        (r, status, kind)
    }
    fn ln(&self) -> (Self, Status) {
        if self.hi <= 0.0 {
            return (Interval::ALL, Status::Certain);
        }
        let s = if self.lo <= 0.0 {
            Status::Maybe
        } else {
            Status::Safe
        };
        let lo = if self.lo <= 0.0 {
            f64::NEG_INFINITY
        } else {
            self.lo.ln()
        };
        (Interval::new(lo, self.hi.ln()), s)
    }
    fn func(op: Op, x: &Self) -> (Self, Status) {
        let x = *x;
        let r = match op {
            Op::Abs => {
                if x.lo >= 0.0 {
                    x
                } else if x.hi <= 0.0 {
                    x.neg()
                } else {
                    Interval::new(0.0, x.lo.abs().max(x.hi))
                }
            }
            Op::Exp => monotone_up(x, f64::exp),
            Op::Floor => monotone_up(x, f64::floor),
            Op::Ceiling => monotone_up(x, f64::ceil),
            Op::Sinh => monotone_up(x, f64::sinh),
            Op::Tanh => monotone_up(x, f64::tanh),
            Op::ArcTan => monotone_up(x, f64::atan),
            Op::ArcSinh => monotone_up(x, f64::asinh),
            Op::Sin | Op::Cos => Interval::new(-1.0, 1.0),
            Op::Cosh => Interval::new(1.0, f64::INFINITY),
            Op::Sech => Interval::new(0.0, 1.0),
            Op::Not => match x.truth() {
                Tri::True => Interval::point(0.0),
                Tri::False => Interval::point(1.0),
                Tri::Unknown => Interval::new(0.0, 1.0),
            },
            // The rest are rare in rate laws: no claim either way.
            _ => Interval::ALL,
        };
        (r, Status::Safe)
    }
    fn min(&self, o: &Self) -> Self {
        Interval::new(self.lo.min(o.lo), self.hi.min(o.hi))
    }
    fn max(&self, o: &Self) -> Self {
        Interval::new(self.lo.max(o.lo), self.hi.max(o.hi))
    }
    fn compare(op: Op, a: &Self, b: &Self) -> Self {
        let t = match op {
            Op::Gt | Op::Lt | Op::Geq | Op::Leq => {
                let (a, b, strict) = match op {
                    Op::Gt => (a, b, true),
                    Op::Lt => (b, a, true),
                    Op::Geq => (a, b, false),
                    _ => (b, a, false),
                };
                // Is a > b (or a >= b)?
                if (strict && a.lo > b.hi) || (!strict && a.lo >= b.hi) {
                    Tri::True
                } else if (strict && a.hi <= b.lo) || (!strict && a.hi < b.lo) {
                    Tri::False
                } else {
                    Tri::Unknown
                }
            }
            Op::Eq | Op::Neq => {
                let eq = if a.is_point() && b.is_point() && a.lo == b.lo {
                    Tri::True
                } else if a.hi < b.lo || b.hi < a.lo {
                    Tri::False
                } else {
                    Tri::Unknown
                };
                match (op, eq) {
                    (Op::Neq, Tri::True) => Tri::False,
                    (Op::Neq, Tri::False) => Tri::True,
                    (_, t) => t,
                }
            }
            _ => Tri::Unknown,
        };
        Self::from_tri(t)
    }
    fn truth(&self) -> Tri {
        if self.lo == 0.0 && self.hi == 0.0 {
            Tri::False
        } else if self.lo > 0.0 || self.hi < 0.0 {
            Tri::True
        } else {
            Tri::Unknown
        }
    }
    fn from_tri(t: Tri) -> Self {
        match t {
            Tri::True => Interval::point(1.0),
            Tri::False => Interval::point(0.0),
            Tri::Unknown => Interval::new(0.0, 1.0),
        }
    }
    fn hull(&self, o: &Self) -> Self {
        Interval::new(self.lo.min(o.lo), self.hi.max(o.hi))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn michaelis_menten_denominator_is_safe() {
        let km = Interval::point(0.5);
        let s = Interval::NON_NEGATIVE;
        let (_, st) = Interval::point(1.0).div(&km.add(&s));
        assert_eq!(st, Status::Safe);
    }

    #[test]
    fn difference_denominator_may_vanish() {
        let s = Interval::NON_NEGATIVE;
        let (_, st) = Interval::point(1.0).div(&s.sub(&s));
        assert_eq!(st, Status::Maybe);
    }

    #[test]
    fn infinite_times_zero() {
        let r = Interval::point(0.0).mul(&Interval::NON_NEGATIVE);
        assert_eq!(r, Interval::point(0.0));
    }

    #[test]
    fn hill_power() {
        let (r, st, _) = Interval::NON_NEGATIVE.pow(&Interval::point(2.5));
        assert_eq!(st, Status::Safe);
        assert_eq!(r.lo, 0.0);
        let (_, st, k) = Interval::NON_NEGATIVE.pow(&Interval::point(-1.0));
        assert_eq!((st, k), (Status::Maybe, HazardKind::ZeroToNegativePower));
    }
}
