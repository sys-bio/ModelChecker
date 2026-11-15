//! A model prepared for analysis: what every id is, what changes over time,
//! and expressions with function calls and assignment rules inlined.

use std::collections::{BTreeMap, HashMap, HashSet};

use sbml_io::{Math, Model, Rule};

use crate::eval::{AVOGADRO, Eval, TIME};
use crate::moiety::Conservation;
use crate::value::Interval;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    Species,
    Compartment,
    Parameter,
    /// In Level 3, a reaction's id stands for its rate.
    Reaction,
    /// A species reference's id stands for its stoichiometry.
    SpeciesReference,
    Function,
}

/// What an analysis may assume about one symbol.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Range {
    /// The values the symbol can take during a simulation.
    pub interval: Interval,
    /// A typical value: its initial value where known.
    pub nominal: f64,
    /// Whether it changes during a simulation.
    pub varies: bool,
}

pub struct View<'m> {
    pub model: &'m Model,
    pub level: u32,
    pub kinds: HashMap<&'m str, Kind>,
    pub functions: HashMap<&'m str, (&'m [String], &'m Math)>,
    pub assignment_rules: HashMap<&'m str, &'m Math>,
    pub rate_rules: HashSet<&'m str>,
    pub event_targets: HashSet<&'m str>,
    pub initial_assignments: HashMap<&'m str, &'m Math>,
    /// Symbols appearing in algebraic rules: they may be set by one.
    pub algebraic: HashSet<String>,
    /// Species changed by at least one reaction (as reactant or product).
    pub in_reactions: HashSet<&'m str>,
    /// Conserved moieties among the species only reactions change.
    pub conservation: Conservation,
    ranges: HashMap<String, Range>,
}

/// The deepest function and rule nesting followed before giving up.
const MAX_DEPTH: usize = 64;

impl<'m> View<'m> {
    pub fn new(model: &'m Model, level: u32) -> Self {
        let mut kinds = HashMap::new();
        for c in &model.compartments {
            kinds.insert(c.id.as_str(), Kind::Compartment);
        }
        for s in &model.species {
            kinds.insert(s.id.as_str(), Kind::Species);
        }
        for p in &model.parameters {
            kinds.insert(p.id.as_str(), Kind::Parameter);
        }
        let mut functions = HashMap::new();
        for f in &model.function_definitions {
            kinds.insert(f.id.as_str(), Kind::Function);
            if let Some(Math::Lambda { args, body }) = f.math.as_ref().map(strip) {
                functions.insert(f.id.as_str(), (args.as_slice(), body.as_ref()));
            }
        }
        let mut in_reactions = HashSet::new();
        for r in &model.reactions {
            if level >= 3 {
                kinds.insert(r.id.as_str(), Kind::Reaction);
            }
            for sr in r.reactants.iter().chain(&r.products) {
                in_reactions.insert(sr.species.as_str());
                if let Some(id) = &sr.id {
                    kinds.insert(id.as_str(), Kind::SpeciesReference);
                }
            }
        }
        let mut assignment_rules = HashMap::new();
        let mut rate_rules = HashSet::new();
        let mut algebraic = HashSet::new();
        for rule in &model.rules {
            match rule {
                Rule::Assignment {
                    variable,
                    math: Some(m),
                    ..
                } => {
                    assignment_rules.insert(variable.as_str(), m);
                }
                Rule::Rate { variable, .. } => {
                    rate_rules.insert(variable.as_str());
                }
                Rule::Algebraic { math: Some(m), .. } => {
                    m.visit_idents(&mut |id| {
                        algebraic.insert(id.to_string());
                    });
                }
                _ => {}
            }
        }
        let event_targets = model
            .events
            .iter()
            .flat_map(|e| e.assignments.iter().map(|a| a.variable.as_str()))
            .collect();
        let initial_assignments = model
            .initial_assignments
            .iter()
            .filter_map(|ia| Some((ia.symbol.as_str(), ia.math.as_ref()?)))
            .collect();
        let mut view = View {
            model,
            level,
            kinds,
            functions,
            assignment_rules,
            rate_rules,
            event_targets,
            initial_assignments,
            algebraic,
            in_reactions,
            conservation: Conservation::new(Vec::new(), Vec::new()),
            ranges: HashMap::new(),
        };
        view.ranges = view.compute_ranges();
        view.conservation = view.compute_conservation();
        view
    }

    fn compute_conservation(&self) -> Conservation {
        let m = self.model;
        let species: Vec<(String, f64)> = m
            .species
            .iter()
            .filter(|s| {
                s.boundary_condition != Some(true)
                    && s.constant != Some(true)
                    && self.in_reactions.contains(s.id.as_str())
                    && !self.changed(&s.id)
                    && !self.assignment_rules.contains_key(s.id.as_str())
            })
            .map(|s| {
                let x0 = self.range(&s.id).map(|r| r.nominal).unwrap_or(f64::NAN);
                (s.id.clone(), if x0.is_finite() { x0 } else { 0.0 })
            })
            .collect();
        let reactions = m
            .reactions
            .iter()
            .map(|r| {
                r.reactants
                    .iter()
                    .map(|sr| (sr.species.clone(), -sr.stoichiometry.unwrap_or(1.0)))
                    .chain(
                        r.products
                            .iter()
                            .map(|sr| (sr.species.clone(), sr.stoichiometry.unwrap_or(1.0))),
                    )
                    .collect()
            })
            .collect();
        Conservation::new(species, reactions)
    }

    /// Whether anything other than an assignment rule changes `id` during
    /// a simulation.
    fn changed(&self, id: &str) -> bool {
        self.rate_rules.contains(id) || self.event_targets.contains(id) || self.algebraic.contains(id)
    }

    fn compute_ranges(&self) -> HashMap<String, Range> {
        let m = self.model;
        let mut r = HashMap::new();
        let fixed = |x: f64| Range {
            interval: Interval::point(x),
            nominal: x,
            varies: false,
        };
        // Pass 1: values as written.
        let mut nominal: BTreeMap<String, f64> = BTreeMap::new();
        for c in &m.compartments {
            nominal.insert(c.id.clone(), c.size.unwrap_or(f64::NAN));
        }
        for s in &m.species {
            let size = m
                .compartments
                .iter()
                .find(|c| c.id == s.compartment)
                .and_then(|c| c.size)
                .unwrap_or(1.0);
            let only_amount = s.has_only_substance_units == Some(true);
            let v = match (s.initial_concentration, s.initial_amount) {
                (_, Some(a)) if only_amount => a,
                (Some(c), _) if only_amount => c * size,
                (Some(c), _) => c,
                (None, Some(a)) => a / size,
                (None, None) => f64::NAN,
            };
            nominal.insert(s.id.clone(), v);
        }
        for p in &m.parameters {
            nominal.insert(p.id.clone(), p.value.unwrap_or(f64::NAN));
        }
        for re in &m.reactions {
            for sr in re.reactants.iter().chain(&re.products) {
                if let Some(id) = &sr.id {
                    nominal.insert(id.clone(), sr.stoichiometry.unwrap_or(1.0));
                }
            }
        }
        // Pass 2: initial assignments and assignment rules, evaluated at
        // the values so far. Twice, so one may use another's result.
        for _ in 0..2 {
            let mut updates = Vec::new();
            let targets = self
                .initial_assignments
                .iter()
                .chain(self.assignment_rules.iter());
            for (id, math) in targets {
                let snapshot = &nominal;
                let env = |k: &str| -> Option<f64> {
                    if k == TIME {
                        return Some(0.0);
                    }
                    if k == AVOGADRO {
                        return Some(6.02214076e23);
                    }
                    snapshot.get(k).copied()
                };
                let expanded = self.expand(math, &HashSet::new());
                let v = Eval::new(&env).eval(&expanded);
                if v.is_finite() {
                    updates.push((id.to_string(), v));
                }
            }
            for (id, v) in updates {
                nominal.insert(id, v);
            }
        }
        // Pass 3: ranges.
        for c in &m.compartments {
            let v = nominal[&c.id];
            let range = if self.changed(&c.id) {
                Range {
                    interval: Interval::NON_NEGATIVE,
                    nominal: if v.is_finite() { v } else { 1.0 },
                    varies: true,
                }
            } else if v.is_finite() {
                fixed(v)
            } else {
                fixed(1.0)
            };
            r.insert(c.id.clone(), range);
        }
        for s in &m.species {
            let v = nominal[&s.id];
            let held = s.boundary_condition == Some(true) || s.constant == Some(true);
            let varies = self.changed(&s.id) || (!held && self.in_reactions.contains(s.id.as_str()));
            let range = if v < 0.0 {
                // A negative initial value (a voltage, a deviation from a
                // reference): this species is not a concentration.
                Range {
                    interval: Interval::ALL,
                    nominal: v,
                    varies,
                }
            } else if varies || !v.is_finite() {
                Range {
                    interval: Interval::NON_NEGATIVE,
                    nominal: if v.is_finite() { v } else { 1.0 },
                    varies: true,
                }
            } else {
                fixed(v)
            };
            r.insert(s.id.clone(), range);
        }
        for p in &m.parameters {
            let v = nominal[&p.id];
            let range = if self.changed(&p.id) || !v.is_finite() {
                let nonneg = !v.is_finite() || v >= 0.0;
                Range {
                    interval: if nonneg {
                        Interval::NON_NEGATIVE
                    } else {
                        Interval::ALL
                    },
                    nominal: if v.is_finite() { v } else { 1.0 },
                    varies: true,
                }
            } else {
                fixed(v)
            };
            r.insert(p.id.clone(), range);
        }
        for (id, v) in nominal.iter() {
            if !r.contains_key(id) {
                r.insert(id.clone(), fixed(*v));
            }
        }
        r.insert(
            TIME.to_string(),
            Range {
                interval: Interval::NON_NEGATIVE,
                nominal: 0.0,
                varies: true,
            },
        );
        r.insert(AVOGADRO.to_string(), fixed(6.02214076e23));
        for re in &m.reactions {
            if self.level >= 3 {
                r.insert(
                    re.id.clone(),
                    Range {
                        interval: Interval::ALL,
                        nominal: f64::NAN,
                        varies: true,
                    },
                );
            }
        }
        r
    }

    /// The range of a global symbol, if it is one.
    pub fn range(&self, id: &str) -> Option<Range> {
        self.ranges.get(id).copied()
    }

    /// `m` with function calls inlined, and assignment-rule variables
    /// replaced by their definitions, except names in `shadow` (a rate
    /// law's local parameters).
    pub fn expand(&self, m: &Math, shadow: &HashSet<&str>) -> Math {
        self.expand_at(m, shadow, true, 0)
    }

    /// `m` with only function calls inlined: what the expression refers to
    /// as written, seen through function definitions.
    pub fn inline_functions(&self, m: &Math) -> Math {
        self.expand_at(m, &HashSet::new(), false, 0)
    }

    fn expand_at(&self, m: &Math, shadow: &HashSet<&str>, rules: bool, depth: usize) -> Math {
        if depth > MAX_DEPTH {
            return Math::Raw("too deeply nested".into());
        }
        let rec = |x: &Math| self.expand_at(x, shadow, rules, depth + 1);
        match m {
            Math::Ident(id) if rules && !shadow.contains(id.as_str()) => {
                match self.assignment_rules.get(id.as_str()) {
                    Some(def) => rec(def),
                    None => m.clone(),
                }
            }
            Math::Call(name, args) => match self.functions.get(name.as_str()) {
                Some((params, body)) if params.len() == args.len() => {
                    let args: Vec<Math> = args.iter().map(&rec).collect();
                    let map: HashMap<&str, &Math> =
                        params.iter().map(String::as_str).zip(args.iter()).collect();
                    rec(&substitute(body, &map))
                }
                _ => Math::Call(name.clone(), args.iter().map(rec).collect()),
            },
            Math::Apply(op, args) => Math::Apply(*op, args.iter().map(rec).collect()),
            Math::CallSymbol(s, args) => Math::CallSymbol(s.clone(), args.iter().map(rec).collect()),
            Math::Root { degree, arg } => Math::Root {
                degree: degree.as_ref().map(|d| Box::new(rec(d))),
                arg: Box::new(rec(arg)),
            },
            Math::Log { base, arg } => Math::Log {
                base: base.as_ref().map(|d| Box::new(rec(d))),
                arg: Box::new(rec(arg)),
            },
            Math::Piecewise { pieces, otherwise } => Math::Piecewise {
                pieces: pieces.iter().map(|(v, c)| (rec(v), rec(c))).collect(),
                otherwise: otherwise.as_ref().map(|o| Box::new(rec(o))),
            },
            Math::Semantics { math, .. } => rec(math),
            _ => m.clone(),
        }
    }
}

fn strip(m: &Math) -> &Math {
    match m {
        Math::Semantics { math, .. } => strip(math),
        _ => m,
    }
}

/// `body` with each `<ci>` named in `map` replaced.
fn substitute(body: &Math, map: &HashMap<&str, &Math>) -> Math {
    let rec = |x: &Math| substitute(x, map);
    match body {
        Math::Ident(id) => match map.get(id.as_str()) {
            Some(v) => (*v).clone(),
            None => body.clone(),
        },
        Math::Apply(op, args) => Math::Apply(*op, args.iter().map(rec).collect()),
        Math::Call(n, args) => Math::Call(n.clone(), args.iter().map(rec).collect()),
        Math::CallSymbol(s, args) => Math::CallSymbol(s.clone(), args.iter().map(rec).collect()),
        Math::Root { degree, arg } => Math::Root {
            degree: degree.as_ref().map(|d| Box::new(rec(d))),
            arg: Box::new(rec(arg)),
        },
        Math::Log { base, arg } => Math::Log {
            base: base.as_ref().map(|d| Box::new(rec(d))),
            arg: Box::new(rec(arg)),
        },
        Math::Piecewise { pieces, otherwise } => Math::Piecewise {
            pieces: pieces.iter().map(|(v, c)| (rec(v), rec(c))).collect(),
            otherwise: otherwise.as_ref().map(|o| Box::new(rec(o))),
        },
        Math::Semantics { math, .. } => rec(math),
        _ => body.clone(),
    }
}

/// Every identifier `m` refers to.
pub fn idents(m: &Math) -> HashSet<String> {
    let mut s = HashSet::new();
    m.visit_idents(&mut |id| {
        s.insert(id.to_string());
    });
    s
}

/// Whether `m` contains MathML that could not be read.
pub fn has_raw(m: &Math) -> bool {
    match m {
        Math::Raw(_) => true,
        Math::Apply(_, a) | Math::Call(_, a) | Math::CallSymbol(_, a) => a.iter().any(has_raw),
        Math::Root { degree: q, arg } | Math::Log { base: q, arg } => {
            q.as_deref().is_some_and(has_raw) || has_raw(arg)
        }
        Math::Piecewise { pieces, otherwise } => {
            pieces.iter().any(|(v, c)| has_raw(v) || has_raw(c))
                || otherwise.as_deref().is_some_and(has_raw)
        }
        Math::Lambda { body, .. } => has_raw(body),
        Math::Semantics { math, .. } => has_raw(math),
        _ => false,
    }
}
