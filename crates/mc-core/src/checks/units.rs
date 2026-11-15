//! Unit checks.
//!
//! UN001 kinetic law units differ from extent/time; UN005 they differ by a
//! compartment's size units (a missing volume factor); UN006 the same
//! suspicion without declared units; UN002/UN003 adding or comparing
//! quantities of different dimensions or scales; UN004 a function argument
//! that is not dimensionless; UN007 a rule or assignment whose units differ
//! from its variable's; UN008 what could not be checked.
//!
//! Units are inferred bottom-up. A symbol without declared units is
//! unknown (a parameter set by an assignment rule takes the rule's units),
//! and an unknown makes a product unknown, so findings arise only where the
//! model says enough. A bare number adapts to what it is added to and is
//! dimensionless in a product.

use std::cell::RefCell;
use std::collections::{BTreeSet, HashMap, HashSet};

use sbml_io::{Math, Op, Rule, Symbol};

use super::scope::{Scope, free_symbols};
use super::{Ctx, locals};
use crate::codes::*;
use crate::finding::{Element, Severity};
use crate::fmt::infix;
use crate::units::{Unit, Units};
use crate::view::{Kind, View};

#[derive(Debug, Clone, Copy)]
enum U {
    Known(Unit),
    /// A number with no units: adapts in sums, dimensionless in products.
    Number,
    Unknown,
}

impl U {
    fn known(self) -> Option<Unit> {
        match self {
            U::Known(u) => Some(u),
            _ => None,
        }
    }
}

struct Issue {
    code: &'static str,
    message: String,
}

struct Infer<'a, 'v, 'm> {
    view: &'v View<'m>,
    units: &'a Units,
    /// Local parameters of the rate law being checked: id to units.
    locals: HashMap<String, Option<Unit>>,
    scope: Scope<'v, 'm>,
    issues: Vec<Issue>,
    /// Inferred units of assignment-rule variables without declared units.
    rule_units: &'a RefCell<HashMap<String, U>>,
    in_progress: &'a RefCell<HashSet<String>>,
}

/// Units a symbol is declared to have (no inference from rules).
fn declared(view: &View, units: &Units, id: &str) -> Option<Unit> {
    let m = view.model;
    match view.kinds.get(id)? {
        Kind::Species => {
            let s = m.species.iter().find(|s| s.id == id)?;
            let substance = match &s.substance_units {
                Some(u) => units.resolve(u)?,
                None => units.substance?,
            };
            if s.has_only_substance_units == Some(true) {
                return Some(substance);
            }
            let c = m.compartments.iter().find(|c| c.id == s.compartment)?;
            if c.spatial_dimensions == Some(0.0) {
                return Some(substance);
            }
            Some(substance.div(&compartment_units(units, c)?))
        }
        Kind::Compartment => {
            let c = m.compartments.iter().find(|c| c.id == id)?;
            compartment_units(units, c)
        }
        Kind::Parameter => {
            let p = m.parameters.iter().find(|p| p.id == id)?;
            units.resolve(p.units.as_deref()?)
        }
        Kind::Reaction => Some(units.extent?.div(&units.time?)),
        Kind::SpeciesReference => Some(Unit::DIMENSIONLESS),
        Kind::Function => None,
    }
}

fn compartment_units(units: &Units, c: &sbml_io::Compartment) -> Option<Unit> {
    if let Some(u) = &c.units {
        return units.resolve(u);
    }
    let d = c.spatial_dimensions.unwrap_or(3.0);
    if d == 3.0 {
        units.volume
    } else if d == 2.0 {
        units.area
    } else if d == 1.0 {
        units.length
    } else if d == 0.0 {
        Some(Unit::DIMENSIONLESS)
    } else {
        None
    }
}

impl Infer<'_, '_, '_> {
    fn issue(&mut self, code: &'static str, message: String) {
        if !self.issues.iter().any(|i| i.message == message) {
            self.issues.push(Issue { code, message });
        }
    }

    fn ident(&mut self, id: &str) -> U {
        if let Some(u) = self.locals.get(id) {
            return u.map_or(U::Unknown, U::Known);
        }
        if let Some(u) = declared(self.view, self.units, id) {
            return U::Known(u);
        }
        // A parameter without units, set by an assignment rule: the rule's.
        if let Some(def) = self.view.assignment_rules.get(id) {
            if let Some(u) = self.rule_units.borrow().get(id) {
                return *u;
            }
            if !self.in_progress.borrow_mut().insert(id.to_string()) {
                return U::Unknown; // a cycle
            }
            let def = self.view.inline_functions(def);
            // Issues inside the rule are reported on the rule itself.
            let saved = std::mem::take(&mut self.issues);
            let saved_locals = std::mem::take(&mut self.locals);
            let u = self.infer(&def);
            self.issues = saved;
            self.locals = saved_locals;
            self.in_progress.borrow_mut().remove(id);
            self.rule_units.borrow_mut().insert(id.to_string(), u);
            return u;
        }
        U::Unknown
    }

    /// Units that must agree: terms of a sum, sides of a comparison,
    /// branches of a piecewise.
    fn agree(&mut self, items: &[(U, &Math)], what: &str) -> U {
        let mut first: Option<(Unit, &Math)> = None;
        let mut any_unknown = false;
        for (u, m) in items {
            match u {
                U::Known(k) => match first {
                    None => first = Some((*k, m)),
                    Some((f, fm)) => {
                        if !f.same_dimensions(k) {
                            let msg = format!(
                                "{what} quantities of different dimensions: `{}` is in {} but `{}` is in {}.",
                                infix(fm),
                                self.units.describe(&f),
                                infix(m),
                                self.units.describe(k)
                            );
                            self.issue(UN_DIMENSION_MISMATCH, msg);
                        } else if !f.same(k) {
                            let msg = format!(
                                "{what} quantities in different units of the same kind: `{}` is in {} but `{}` is in {} (a factor of {}).",
                                infix(fm),
                                self.units.describe(&f),
                                infix(m),
                                self.units.describe(k),
                                crate::sample::fmt_num(f.factor(k))
                            );
                            self.issue(UN_SCALE_MISMATCH, msg);
                        }
                    }
                },
                U::Unknown => any_unknown = true,
                U::Number => {}
            }
        }
        match first {
            Some((u, _)) => U::Known(u),
            None if any_unknown => U::Unknown,
            None => U::Number,
        }
    }

    fn dimensionless_arg(&mut self, u: U, arg: &Math, what: &str) {
        if let U::Known(k) = u
            && !k.is_dimensionless()
        {
            let msg = format!(
                "The argument of {what} should be dimensionless, but `{}` is in {}.",
                infix(arg),
                self.units.describe(&k)
            );
            self.issue(UN_FUNCTION_ARGUMENT, msg);
        }
    }

    /// The numeric value of an exponent, at the model's values.
    fn value(&self, m: &Math) -> Option<f64> {
        let (v, _) = self.scope.at(m, &[], &[]);
        v.is_finite().then_some(v)
    }

    fn infer(&mut self, m: &Math) -> U {
        match m {
            Math::Number { units, .. } => match units {
                Some(u) => self.units.resolve(u).map_or(U::Unknown, U::Known),
                None => U::Number,
            },
            Math::Ident(id) => self.ident(id),
            Math::Symbol(s) => {
                if s.url == Symbol::TIME {
                    self.units.time.map_or(U::Unknown, U::Known)
                } else if s.url == Symbol::AVOGADRO {
                    U::Known(crate::units::base_unit("mole").unwrap().pow(-1.0))
                } else {
                    U::Unknown
                }
            }
            Math::Constant(_) => U::Number,
            Math::Apply(op, args) => self.apply(*op, args),
            Math::Call(..) => U::Unknown,
            Math::CallSymbol(s, args) => {
                if s.url == Symbol::DELAY && !args.is_empty() {
                    for a in &args[1..] {
                        self.infer(a);
                    }
                    self.infer(&args[0])
                } else if s.url == Symbol::RATE_OF && args.len() == 1 {
                    match (self.infer(&args[0]), self.units.time) {
                        (U::Known(u), Some(t)) => U::Known(u.div(&t)),
                        _ => U::Unknown,
                    }
                } else {
                    U::Unknown
                }
            }
            Math::Root { degree, arg } => {
                let u = self.infer(arg);
                let d = match degree {
                    Some(d) => self.value(d),
                    None => Some(2.0),
                };
                match (u, d) {
                    (U::Known(k), _) if k.is_dimensionless() => U::Known(k),
                    (U::Known(k), Some(d)) if d != 0.0 => U::Known(k.pow(1.0 / d)),
                    (U::Number, _) => U::Number,
                    _ => U::Unknown,
                }
            }
            Math::Log { base, arg } => {
                if let Some(b) = base {
                    let bu = self.infer(b);
                    self.dimensionless_arg(bu, b, "a logarithm's base");
                }
                let u = self.infer(arg);
                self.dimensionless_arg(u, arg, "log");
                U::Known(Unit::DIMENSIONLESS)
            }
            Math::Piecewise { pieces, otherwise } => {
                let mut values: Vec<(U, &Math)> = Vec::new();
                for (v, c) in pieces {
                    self.infer(c);
                    values.push((self.infer(v), v));
                }
                if let Some(o) = otherwise {
                    values.push((self.infer(o), o));
                }
                self.agree(&values, "A piecewise expression has branches with")
            }
            Math::Lambda { body, .. } => self.infer(body),
            Math::Semantics { math, .. } => self.infer(math),
            Math::Raw(_) => U::Unknown,
        }
    }

    fn apply(&mut self, op: Op, args: &[Math]) -> U {
        match op {
            Op::Plus | Op::Minus | Op::Min | Op::Max => {
                let us: Vec<(U, &Math)> = args.iter().map(|a| (self.infer(a), a)).collect();
                let what = match op {
                    Op::Plus => "Adds",
                    Op::Minus => "Subtracts",
                    _ => "Compares",
                };
                self.agree(&us, what)
            }
            Op::Times | Op::Divide => {
                let us: Vec<U> = args.iter().map(|a| self.infer(a)).collect();
                if us.iter().any(|u| matches!(u, U::Unknown)) {
                    return U::Unknown;
                }
                // A number with no units in a product may carry units of its
                // own (a rate constant written as a number, a conversion
                // factor). Only small whole numbers and 1/2, which are
                // almost always stoichiometric, count as dimensionless.
                let plain = |a: &Math| match a {
                    Math::Number { value, units: None } => {
                        let v = value.value().abs();
                        (v.fract() == 0.0 && v <= 1000.0) || v == 0.5
                    }
                    _ => true,
                };
                if !args.iter().all(plain) {
                    return U::Unknown;
                }
                if us.iter().all(|u| matches!(u, U::Number)) {
                    return U::Number;
                }
                let mut acc = Unit::DIMENSIONLESS;
                for (i, u) in us.iter().enumerate() {
                    let k = u.known().unwrap_or(Unit::DIMENSIONLESS);
                    acc = if op == Op::Divide && i > 0 {
                        acc.div(&k)
                    } else {
                        acc.mul(&k)
                    };
                }
                U::Known(acc)
            }
            Op::Power => {
                let [a, b] = args else { return U::Unknown };
                let base = self.infer(a);
                let e = self.infer(b);
                self.dimensionless_arg(e, b, "an exponent");
                match base {
                    U::Known(k) if k.is_dimensionless() => U::Known(k),
                    U::Known(k) => match self.value(b) {
                        Some(p) => U::Known(k.pow(p)),
                        None => U::Unknown,
                    },
                    other => other,
                }
            }
            Op::Exp | Op::Ln | Op::Factorial | Op::Sin | Op::Cos | Op::Tan | Op::Sec
            | Op::Csc | Op::Cot | Op::Sinh | Op::Cosh | Op::Tanh | Op::Sech | Op::Csch
            | Op::Coth | Op::ArcSin | Op::ArcCos | Op::ArcTan | Op::ArcSec | Op::ArcCsc
            | Op::ArcCot | Op::ArcSinh | Op::ArcCosh | Op::ArcTanh | Op::ArcSech
            | Op::ArcCsch | Op::ArcCoth => {
                for a in args {
                    let u = self.infer(a);
                    self.dimensionless_arg(u, a, op.tag());
                }
                U::Known(Unit::DIMENSIONLESS)
            }
            Op::Abs | Op::Floor | Op::Ceiling => match args {
                [a] => self.infer(a),
                _ => U::Unknown,
            },
            Op::Eq | Op::Neq | Op::Gt | Op::Lt | Op::Geq | Op::Leq => {
                let us: Vec<(U, &Math)> = args.iter().map(|a| (self.infer(a), a)).collect();
                self.agree(&us, "Compares");
                U::Known(Unit::DIMENSIONLESS)
            }
            Op::And | Op::Or | Op::Xor | Op::Not | Op::Implies => {
                for a in args {
                    self.infer(a);
                }
                U::Known(Unit::DIMENSIONLESS)
            }
            Op::Quotient | Op::Rem => {
                for a in args {
                    self.infer(a);
                }
                U::Unknown
            }
        }
    }
}

pub fn check(ctx: &mut Ctx) {
    let view = ctx.view;
    let m = view.model;
    let units = Units::new(m, view.level);
    let rule_units = RefCell::new(HashMap::new());
    let in_progress = RefCell::new(HashSet::new());
    let declared_any = units.any_declared();

    let new_infer = |locals: HashMap<String, Option<Unit>>, values: HashMap<String, f64>| {
        let mut scope = Scope::new(view);
        scope.locals = values;
        Infer {
            view,
            units: &units,
            locals,
            scope,
            issues: Vec::new(),
            rule_units: &rule_units,
            in_progress: &in_progress,
        }
    };

    let expected_rate = match (units.extent, units.time) {
        (Some(e), Some(t)) => Some(e.div(&t)),
        _ => None,
    };
    let mut undetermined = 0usize;
    let mut laws = 0usize;
    let mut missing_units: BTreeSet<String> = BTreeSet::new();
    // Per-reaction rate findings that may share a cause: (index of the
    // finding, its group, the reaction).
    let mut grouped: Vec<(usize, Group, String)> = Vec::new();

    for r in &m.reactions {
        let Some(math) = r.kinetic_law.as_ref().and_then(|k| k.math.as_ref()) else {
            continue;
        };
        laws += 1;
        let el = Element::Reaction(r.id.clone());
        let local_units: HashMap<String, Option<Unit>> = r
            .kinetic_law
            .iter()
            .flat_map(|k| &k.parameters)
            .map(|p| (p.id.clone(), p.units.as_deref().and_then(|u| units.resolve(u))))
            .collect();
        let inlined = view.inline_functions(math);
        let mut inf = new_infer(local_units.clone(), locals(r));
        let u = inf.infer(&inlined);
        let issues = std::mem::take(&mut inf.issues);
        for i in issues {
            ctx.add(i.code, Severity::Warning, el.clone(), with_law(&i.message, &inlined), None);
        }
        let mut volume_reported = false;
        match (u, expected_rate) {
            (U::Known(got), Some(want)) if !got.same(&want) => {
                let g = check_rate_units(ctx, &units, r, &inlined, got, want);
                volume_reported = g.code == UN_MISSING_VOLUME;
                grouped.push((ctx.findings.len() - 1, g, r.id.clone()));
            }
            (U::Known(_), _) => {}
            _ => {
                undetermined += 1;
                for id in free_symbols(&inlined) {
                    let unitless = match local_units.get(&id) {
                        Some(u) => u.is_none(),
                        None => {
                            view.kinds.get(id.as_str()) == Some(&Kind::Parameter)
                                && declared(view, &units, &id).is_none()
                                && !view.assignment_rules.contains_key(id.as_str())
                        }
                    };
                    if unitless {
                        missing_units.insert(id);
                    }
                }
            }
        }
        if !volume_reported
            && !matches!(u, U::Known(_))
            && let Some(g) = volume_heuristic(ctx, r, math)
        {
            grouped.push((ctx.findings.len() - 1, g, r.id.clone()));
        }
    }
    merge_groups(ctx, grouped);

    // Rules and assignments against their variables' declared units.
    let check_target = |ctx: &mut Ctx, el: Element, var: &str, math: &Math, per_time: bool, what: &str| {
        let Some(mut want) = declared(view, &units, var) else {
            return;
        };
        if per_time {
            let Some(t) = units.time else { return };
            want = want.div(&t);
        }
        let inlined = view.inline_functions(math);
        let mut inf = new_infer(HashMap::new(), HashMap::new());
        let got = inf.infer(&inlined);
        for i in std::mem::take(&mut inf.issues) {
            ctx.add(i.code, Severity::Warning, el.clone(), with_law(&i.message, &inlined), None);
        }
        if let U::Known(g) = got
            && !g.same(&want)
        {
            ctx.add(
                UN_ASSIGNMENT,
                Severity::Warning,
                el,
                format!(
                    "The {what} for `{var}` has units {}, but {} {}.\nExpression: {}",
                    units.describe(&g),
                    if per_time { "the variable's units per time are" } else { "the variable is in" },
                    units.describe(&want),
                    infix(&inlined)
                ),
                None,
            );
        }
    };
    for rule in &m.rules {
        match rule {
            Rule::Assignment { variable, math: Some(math), .. } => {
                check_target(ctx, Element::Rule(variable.clone()), variable, math, false, "assignment rule");
            }
            Rule::Rate { variable, math: Some(math), .. } => {
                check_target(ctx, Element::Rule(variable.clone()), variable, math, true, "rate rule");
            }
            _ => {}
        }
    }
    for ia in &m.initial_assignments {
        if let Some(math) = &ia.math {
            check_target(ctx, Element::InitialAssignment(ia.symbol.clone()), &ia.symbol, math, false, "initial assignment");
        }
    }
    for (i, e) in m.events.iter().enumerate() {
        let id = e.id.clone().unwrap_or_else(|| format!("#{}", i + 1));
        for a in &e.assignments {
            if let Some(math) = &a.math {
                check_target(ctx, Element::Event(id.clone()), &a.variable, math, false, "event assignment");
            }
        }
    }

    // What could not be checked.
    if laws > 0 {
        if !declared_any {
            ctx.add(
                UN_COVERAGE,
                Severity::Info,
                Element::Model,
                "The model declares no units, so the units of its rate laws could not be checked.",
                None,
            );
        } else if expected_rate.is_none() {
            ctx.add(
                UN_COVERAGE,
                Severity::Info,
                Element::Model,
                "The model does not declare its extent and time units, so the units of its rate laws could not be checked against them.",
                None,
            );
        } else if undetermined > 0 {
            let mut names: Vec<&String> = missing_units.iter().collect();
            let more = names.len().saturating_sub(10);
            names.truncate(10);
            let list = names.iter().map(|s| s.as_str()).collect::<Vec<_>>().join(", ");
            ctx.add(
                UN_COVERAGE,
                Severity::Info,
                Element::Model,
                format!(
                    "The units of {undetermined} of {laws} rate laws could not be determined{}.",
                    if list.is_empty() {
                        String::new()
                    } else if more > 0 {
                        format!(", because these have no units: {list} and {more} more")
                    } else {
                        format!(", because these have no units: {list}")
                    }
                ),
                None,
            );
        }
    }
}

fn with_law(message: &str, m: &Math) -> String {
    format!("{message}\nExpression: {}", infix(m))
}

/// What a per-reaction rate finding has in common with others of the same
/// cause, for merging them.
struct Group {
    code: &'static str,
    key: String,
    /// "N rate laws ..." continues with this.
    summary: String,
    /// Why they may share a cause.
    hint: String,
}

/// The fewest findings with one cause that are merged into one.
const MERGE_AT: usize = 3;

/// Replaces each group of at least [`MERGE_AT`] findings with one cause by
/// a single model-level finding that lists the reactions.
fn merge_groups(ctx: &mut Ctx, grouped: Vec<(usize, Group, String)>) {
    let mut by_key: Vec<(&Group, Vec<(usize, &str)>)> = Vec::new();
    for (i, g, r) in &grouped {
        match by_key
            .iter_mut()
            .find(|(k, _)| k.code == g.code && k.key == g.key)
        {
            Some((_, v)) => v.push((*i, r)),
            None => by_key.push((g, vec![(*i, r)])),
        }
    }
    let mut remove: Vec<usize> = Vec::new();
    let mut merged = Vec::new();
    for (g, members) in by_key {
        if members.len() < MERGE_AT {
            continue;
        }
        remove.extend(members.iter().map(|(i, _)| *i));
        let mut ids: Vec<&str> = members.iter().map(|(_, r)| *r).collect();
        let n = ids.len();
        let more = n.saturating_sub(12);
        ids.truncate(12);
        let list = if more > 0 {
            format!("{} and {more} more", ids.join(", "))
        } else {
            ids.join(", ")
        };
        merged.push(crate::finding::Finding {
            code: g.code,
            severity: Severity::Warning,
            element: Element::Model,
            message: format!("{n} rate laws {}: {list}. {}", g.summary, g.hint),
            witness: None,
            reaction: None,
            line: None,
        });
    }
    remove.sort_unstable();
    for i in remove.into_iter().rev() {
        ctx.findings.remove(i);
    }
    ctx.findings.extend(merged);
}

/// UN001 or UN005 for a rate law whose units are known and wrong.
fn check_rate_units(
    ctx: &mut Ctx,
    units: &Units,
    r: &sbml_io::Reaction,
    law: &Math,
    got: Unit,
    want: Unit,
) -> Group {
    let view = ctx.view;
    let el = Element::Reaction(r.id.clone());
    // Is it short of exactly one compartment's size units?
    let comps: BTreeSet<&str> = r
        .reactants
        .iter()
        .chain(&r.products)
        .filter_map(|sr| view.model.species.iter().find(|s| s.id == sr.species))
        .map(|s| s.compartment.as_str())
        .collect();
    for c in comps {
        let Some(cu) = declared(view, units, c) else { continue };
        if !cu.is_dimensionless() && got.mul(&cu).same(&want) {
            ctx.add(
                UN_MISSING_VOLUME,
                Severity::Warning,
                el,
                format!(
                    "The rate law is in {} (concentration per time), but a kinetic law must give \
                     {} (amount per time). It is probably missing a factor of the compartment \
                     size, `{c}`.\nRate law: {}",
                    units.describe(&got),
                    units.describe(&want),
                    infix(law)
                ),
                None,
            );
            return Group {
                code: UN_MISSING_VOLUME,
                key: c.to_string(),
                summary: format!(
                    "are in concentration per time, short of amount per time by the size of compartment `{c}`"
                ),
                hint: format!("Multiply each by `{c}`."),
            };
        }
    }
    let note = if got.same_dimensions(&want) {
        format!(
            " The dimensions agree but the scale differs by a factor of {}.",
            crate::sample::fmt_num(got.factor(&want))
        )
    } else {
        String::new()
    };
    ctx.add(
        UN_RATE_UNITS,
        Severity::Warning,
        el,
        format!(
            "The rate law is in {}, but a kinetic law must be in extent per time, {}.{note}\nRate law: {}",
            units.describe(&got),
            units.describe(&want),
            infix(law)
        ),
        None,
    );
    // Rate laws off by the same ratio probably share one cause.
    let ratio = got.div(&want);
    let time_only = (0..8).all(|i| i == 2 || ratio.exp[i].abs() < 1e-9) && ratio.exp[2].abs() > 1e-9;
    let hint = if ratio.is_dimensionless() {
        format!(
            "They are all off by the same factor, {}, which points to a unit definition (often the substance or time unit) being off by that factor, not to the rate laws themselves.",
            crate::sample::fmt_num(got.factor(&want))
        )
    } else if time_only {
        "They all differ by a time dimension, which points to rate constants declared with the wrong time units, or with none."
            .to_string()
    } else {
        "They all differ in the same way, so the cause is probably shared: a unit definition, or a parameter they all use."
            .to_string()
    };
    Group {
        code: UN_RATE_UNITS,
        key: ratio.to_string(),
        summary: if ratio.is_dimensionless() {
            format!(
                "are off from extent per time ({}) by the same factor, {}",
                units.describe(&want),
                crate::sample::fmt_num(got.factor(&want))
            )
        } else {
            format!(
                "are all in {} where extent per time ({}) is needed",
                units.describe(&got),
                units.describe(&want)
            )
        },
        hint,
    }
}

/// UN006: without units to go on, a rate law that uses concentrations of
/// species in a compartment whose size is not 1, and never mentions a
/// compartment, is probably missing the volume factor.
fn volume_heuristic(ctx: &mut Ctx, r: &sbml_io::Reaction, math: &Math) -> Option<Group> {
    let view = ctx.view;
    let m = view.model;
    let shadow = super::shadow(&locals(r)).iter().map(|s| s.to_string()).collect::<Vec<_>>();
    let shadow: HashSet<&str> = shadow.iter().map(String::as_str).collect();
    let expanded = view.expand(math, &shadow);
    let deps: HashSet<String> = free_symbols(&expanded).into_iter().collect();
    if deps
        .iter()
        .any(|d| view.kinds.get(d.as_str()) == Some(&Kind::Compartment))
    {
        return None;
    }
    // Species the rate depends on, as concentrations, in a compartment of
    // size other than 1.
    let mut culprit: Option<(&str, &str, f64)> = None;
    for sr in r.reactants.iter().chain(&r.products) {
        let Some(s) = m.species.iter().find(|s| s.id == sr.species) else { continue };
        if !deps.contains(&s.id) || s.has_only_substance_units == Some(true) {
            continue;
        }
        let Some(c) = m.compartments.iter().find(|c| c.id == s.compartment) else { continue };
        if c.spatial_dimensions == Some(0.0) {
            continue;
        }
        let size = view.range(&c.id).map(|r| r.nominal).unwrap_or(1.0);
        if size.is_finite() && (size - 1.0).abs() > 1e-9 {
            culprit = Some((s.id.as_str(), c.id.as_str(), size));
            break;
        }
    }
    let (s, c, size) = culprit?;
    ctx.add(
        UN_VOLUME_HEURISTIC,
        Severity::Warning,
        Element::Reaction(r.id.clone()),
        format!(
            "The rate law uses the concentration of {s}, in compartment `{c}` of size {}, but \
             never refers to a compartment. A kinetic law gives amount per time, so it usually \
             needs a factor of the compartment size; without it the rate is off by a factor of \
             {}.\nRate law: {}",
            crate::sample::fmt_num(size),
            crate::sample::fmt_num(size),
            infix(&expanded)
        ),
        None,
    );
    Some(Group {
        code: UN_VOLUME_HEURISTIC,
        key: c.to_string(),
        summary: format!(
            "use concentrations of species in compartment `{c}` (size {}) but never refer to a compartment",
            crate::sample::fmt_num(size)
        ),
        hint: format!(
            "Kinetic laws give amount per time, so each usually needs a factor of `{c}`, unless its rate constants already include the volume."
        ),
    })
}
