//! Rules and events.
//!
//! - RU001 assignment rules depending on each other in a loop (through
//!   reaction rates too, in Level 3); RU002 the same among initial
//!   assignments;
//! - RU003 a variable set more than once: two rules, a rule and an initial
//!   assignment, a rule and an event, or a species both changed by
//!   reactions and set by a rule;
//! - RU004 something marked constant that a rule or event sets;
//! - RU005 algebraic rules (many simulators do not support them);
//! - RU006 an event trigger already true at the start;
//! - RU007 an event whose trigger can never become true;
//! - RU008 an event with no assignments; RU009 a negative delay.

use std::collections::{BTreeMap, BTreeSet, HashSet};

use sbml_io::{Event, Math, Rule};

use super::Ctx;
use super::scope::{Scope, free_symbols};
use crate::codes::*;
use crate::finding::{Element, Severity};
use crate::fmt::infix;
use crate::sample::fmt_num;
use crate::view::Kind;

pub fn check(ctx: &mut Ctx) {
    loops(ctx);
    set_more_than_once(ctx);
    algebraic(ctx);
    events(ctx);
}

// ------------------------------------------------------------- loops

/// Strongly connected components with a cycle (more than one node, or a
/// node depending on itself), each as its nodes in order.
fn cycles(graph: &BTreeMap<String, BTreeSet<String>>) -> Vec<Vec<String>> {
    // Tarjan's algorithm, iteratively.
    struct T<'g> {
        g: &'g BTreeMap<String, BTreeSet<String>>,
        index: BTreeMap<&'g str, usize>,
        low: BTreeMap<&'g str, usize>,
        on: HashSet<&'g str>,
        stack: Vec<&'g str>,
        next: usize,
        out: Vec<Vec<String>>,
    }
    fn visit<'g>(t: &mut T<'g>, v: &'g str) {
        t.index.insert(v, t.next);
        t.low.insert(v, t.next);
        t.next += 1;
        t.stack.push(v);
        t.on.insert(v);
        if let Some(ws) = t.g.get(v) {
            for w in ws {
                let w = w.as_str();
                if !t.g.contains_key(w) {
                    continue;
                }
                if !t.index.contains_key(w) {
                    visit(t, w);
                    let lw = t.low[w];
                    let lv = t.low.get_mut(v).unwrap();
                    *lv = (*lv).min(lw);
                } else if t.on.contains(w) {
                    let iw = t.index[w];
                    let lv = t.low.get_mut(v).unwrap();
                    *lv = (*lv).min(iw);
                }
            }
        }
        if t.low[v] == t.index[v] {
            let mut comp = Vec::new();
            while let Some(w) = t.stack.pop() {
                t.on.remove(w);
                comp.push(w.to_string());
                if w == v {
                    break;
                }
            }
            let self_loop = t.g.get(v).is_some_and(|ws| ws.contains(v));
            if comp.len() > 1 || self_loop {
                comp.reverse();
                t.out.push(comp);
            }
        }
    }
    let mut t = T {
        g: graph,
        index: BTreeMap::new(),
        low: BTreeMap::new(),
        on: HashSet::new(),
        stack: Vec::new(),
        next: 0,
        out: Vec::new(),
    };
    for v in graph.keys() {
        if !t.index.contains_key(v.as_str()) {
            visit(&mut t, v);
        }
    }
    t.out
}

fn describe_cycle(c: &[String]) -> String {
    if c.len() == 1 {
        return format!("{} depends on itself", c[0]);
    }
    let mut parts: Vec<String> = c.windows(2).map(|w| format!("{} depends on {}", w[0], w[1])).collect();
    parts.push(format!("which depends on {}", c[0]));
    parts.join(", ")
}

/// RU001 and RU002.
fn loops(ctx: &mut Ctx) {
    let view = ctx.view;
    let m = view.model;
    let deps = |math: &Math| -> BTreeSet<String> {
        free_symbols(&view.inline_functions(math)).into_iter().collect()
    };
    // Assignment rules, and in Level 3 reactions (whose ids are their rates).
    let mut graph: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for rule in &m.rules {
        if let Rule::Assignment { variable, math: Some(x), .. } = rule {
            graph.insert(variable.clone(), deps(x));
        }
    }
    if view.level >= 3 {
        for r in &m.reactions {
            if let Some(x) = r.kinetic_law.as_ref().and_then(|k| k.math.as_ref()) {
                let locals: HashSet<&str> = r
                    .kinetic_law
                    .iter()
                    .flat_map(|k| &k.parameters)
                    .map(|p| p.id.as_str())
                    .collect();
                let d = deps(x).into_iter().filter(|s| !locals.contains(s.as_str())).collect();
                graph.insert(r.id.clone(), d);
            }
        }
    }
    for c in cycles(&graph) {
        // A loop only among reactions cannot happen; one must hold a rule.
        let el = c
            .iter()
            .find(|v| view.assignment_rules.contains_key(v.as_str()))
            .cloned()
            .unwrap_or_else(|| c[0].clone());
        ctx.add(
            RU_RULE_LOOP,
            Severity::Error,
            Element::Rule(el),
            format!(
                "Assignment rules depend on each other in a loop: {}. Each would need the \
                 other's value first, so none can be computed; SBML does not allow this.",
                describe_cycle(&c)
            ),
            None,
        );
    }

    // Initial assignments, which are worked out together with assignment
    // rules at the start.
    let mut init: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for ia in &m.initial_assignments {
        if let Some(x) = &ia.math {
            init.insert(ia.symbol.clone(), deps(x));
        }
    }
    if init.is_empty() {
        return;
    }
    for rule in &m.rules {
        if let Rule::Assignment { variable, math: Some(x), .. } = rule {
            init.entry(variable.clone()).or_insert_with(|| deps(x));
        }
    }
    for c in cycles(&init) {
        // Loops among rules alone are reported above.
        let Some(ia) = c.iter().find(|v| view.initial_assignments.contains_key(v.as_str())) else {
            continue;
        };
        ctx.add(
            RU_INITIAL_LOOP,
            Severity::Error,
            Element::InitialAssignment(ia.clone()),
            format!(
                "Initial values depend on each other in a loop: {}. None of them can be \
                 worked out; SBML does not allow this.",
                describe_cycle(&c)
            ),
            None,
        );
    }
}

// ----------------------------------------------------- set more than once

/// RU003 and RU004.
fn set_more_than_once(ctx: &mut Ctx) {
    let view = ctx.view;
    let m = view.model;
    let mut rules: BTreeMap<&str, Vec<&'static str>> = BTreeMap::new();
    for rule in &m.rules {
        match rule {
            Rule::Assignment { variable, .. } => rules.entry(variable).or_default().push("an assignment rule"),
            Rule::Rate { variable, .. } => rules.entry(variable).or_default().push("a rate rule"),
            Rule::Algebraic { .. } => {}
        }
    }
    let element = |id: &str| -> Element {
        match view.kinds.get(id) {
            Some(Kind::Species) => Element::Species(id.to_string()),
            Some(Kind::Compartment) => Element::Compartment(id.to_string()),
            _ => Element::Parameter(id.to_string()),
        }
    };
    for (v, kinds) in &rules {
        if kinds.len() > 1 {
            ctx.add(
                RU_SET_TWICE,
                Severity::Error,
                element(v),
                format!("`{v}` is set by {} rules ({}); a variable may have only one.", kinds.len(), kinds.join(" and ")),
                None,
            );
        }
    }
    let assigned: HashSet<&str> = m
        .rules
        .iter()
        .filter_map(|r| match r {
            Rule::Assignment { variable, .. } => Some(variable.as_str()),
            _ => None,
        })
        .collect();
    for ia in &m.initial_assignments {
        if assigned.contains(ia.symbol.as_str()) {
            ctx.add(
                RU_SET_TWICE,
                Severity::Error,
                element(&ia.symbol),
                format!(
                    "`{}` has both an assignment rule and an initial assignment. The rule sets \
                     it at all times, the start included, so the initial assignment conflicts; \
                     SBML does not allow both.",
                    ia.symbol
                ),
                None,
            );
        }
    }
    let mut reported_event: HashSet<&str> = HashSet::new();
    for (i, e) in m.events.iter().enumerate() {
        for a in &e.assignments {
            if assigned.contains(a.variable.as_str()) && reported_event.insert(a.variable.as_str()) {
                ctx.add(
                    RU_SET_TWICE,
                    Severity::Error,
                    event_element(e, i),
                    format!(
                        "The event sets `{}`, which an assignment rule also sets at all times, so \
                         the event's change would be overwritten at once; SBML does not allow \
                         this.",
                        a.variable
                    ),
                    None,
                );
            }
        }
    }
    // A species changed by reactions and also set by a rule.
    for s in &m.species {
        if s.boundary_condition == Some(true) {
            continue;
        }
        let Some(kinds) = rules.get(s.id.as_str()) else { continue };
        if !view.in_reactions.contains(s.id.as_str()) {
            continue;
        }
        ctx.add(
            RU_SET_TWICE,
            Severity::Error,
            Element::Species(s.id.clone()),
            format!(
                "{} is changed by reactions and also set by {}. A species set by a rule must be \
                 a boundary species if it takes part in reactions; otherwise the reactions and \
                 the rule both claim it.",
                s.id,
                kinds[0]
            ),
            None,
        );
    }

    // RU004: constants that change.
    let constant = |id: &str| -> bool {
        m.parameters.iter().any(|p| p.id == id && p.constant == Some(true))
            || m.compartments.iter().any(|c| c.id == id && c.constant == Some(true))
            || m.species.iter().any(|s| s.id == id && s.constant == Some(true))
    };
    for (v, kinds) in &rules {
        if constant(v) {
            ctx.add(
                RU_CONSTANT_CHANGED,
                Severity::Error,
                element(v),
                format!("`{v}` is marked constant, but {} sets it.", kinds[0]),
                None,
            );
        }
    }
    let mut seen: HashSet<&str> = HashSet::new();
    for (i, e) in m.events.iter().enumerate() {
        for a in &e.assignments {
            if constant(&a.variable) && seen.insert(a.variable.as_str()) {
                ctx.add(
                    RU_CONSTANT_CHANGED,
                    Severity::Error,
                    event_element(e, i),
                    format!("The event changes `{}`, which is marked constant.", a.variable),
                    None,
                );
            }
        }
    }
}

fn event_element(e: &Event, i: usize) -> Element {
    Element::Event(e.id.clone().unwrap_or_else(|| format!("#{}", i + 1)))
}

// --------------------------------------------------------- algebraic

/// RU005.
fn algebraic(ctx: &mut Ctx) {
    let n = ctx
        .view
        .model
        .rules
        .iter()
        .filter(|r| matches!(r, Rule::Algebraic { .. }))
        .count();
    if n > 0 {
        ctx.add(
            RU_ALGEBRAIC,
            Severity::Info,
            Element::Model,
            format!(
                "The model has {n} algebraic rule{} (a constraint 0 = f(...) the simulator must \
                 solve for). Many simulators do not support algebraic rules.",
                if n == 1 { "" } else { "s" }
            ),
            None,
        );
    }
}

// ------------------------------------------------------------ events

/// RU006 to RU009.
fn events(ctx: &mut Ctx) {
    let view = ctx.view;
    let m = view.model;
    let scope = Scope::new(view);
    for (i, e) in m.events.iter().enumerate() {
        let el = event_element(e, i);
        if e.assignments.is_empty() {
            ctx.add(
                RU_EVENT_EMPTY,
                Severity::Info,
                el.clone(),
                "The event changes nothing: it has no event assignments.",
                None,
            );
        }
        if let Some(d) = e.delay.as_ref().and_then(|d| d.math.as_ref()) {
            let d = view.expand(d, &HashSet::new());
            let (v, _) = scope.at(&d, &[], &[]);
            if v.is_finite() && v < 0.0 {
                ctx.add(
                    RU_NEGATIVE_DELAY,
                    Severity::Warning,
                    el.clone(),
                    format!(
                        "The event's delay is {} at the start; a delay must not be negative.\nExpression: {}",
                        fmt_num(v),
                        infix(&d)
                    ),
                    None,
                );
            }
        }
        let Some(trigger) = e.trigger.as_ref().and_then(|t| t.math.as_ref()) else {
            continue;
        };
        let trigger = view.expand(trigger, &HashSet::new());
        let (v, _) = scope.at(&trigger, &[], &[]);
        if !v.is_finite() {
            continue;
        }
        let true_at_start = v != 0.0;
        // Level 2 behaves as if the trigger were true just before the start.
        let initial_value = e.trigger.as_ref().and_then(|t| t.initial_value).unwrap_or(true);
        let changes = free_symbols(&trigger)
            .iter()
            .any(|s| scope.range(s).is_some_and(|r| r.varies));
        if !changes {
            // Nothing in the trigger ever changes: it is true or false for good.
            if !true_at_start || initial_value {
                ctx.add(
                    RU_EVENT_NEVER,
                    Severity::Warning,
                    el,
                    format!(
                        "The event can never fire: its trigger depends on nothing that changes \
                         during a simulation, and it {}.\nTrigger: {}",
                        if true_at_start {
                            "is already true at the start (an event fires only when its trigger becomes true)"
                        } else {
                            "is false"
                        },
                        infix(&trigger)
                    ),
                    None,
                );
                continue;
            }
        }
        if true_at_start {
            let consequence = if initial_value {
                "Because the trigger counts as already true before the start (initialValue is \
                 true, or this is Level 2), the event does not fire at time 0; it fires only \
                 after the trigger becomes false and then true again."
            } else {
                "Because initialValue is false, the event fires at once, at time 0."
            };
            ctx.add(
                RU_TRIGGER_TRUE_AT_START,
                Severity::Warning,
                el,
                format!(
                    "The event's trigger is already true at the initial values. {consequence} \
                     Check that this is what is meant.\nTrigger: {}",
                    infix(&trigger)
                ),
                None,
            );
        }
    }
}
