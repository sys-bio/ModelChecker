//! Applying an edit of the Antimony text to the original SBML model,
//! without losing anything.
//!
//! The original document stays the master copy. The Antimony text is a
//! view of it, so an edit is applied as a three-way merge:
//!
//! - `base`: the model read back from the text as first written,
//! - `mine`: the model read back from the edited text,
//! - `orig`: the original model.
//!
//! For every element (by id), a field is changed in `orig` only if it
//! differs between `base` and `mine`, which is exactly when the user
//! changed it. Everything else on the element, and everything the text
//! cannot express (annotations, notes, names, model-wide units, rate-law
//! SBO terms, constraints, layout), is kept as it was. Units and SBO terms
//! the text does show are edited like any other field. An element missing from
//! the edited text was deleted by the user; one missing from `base` was
//! added. An unedited text therefore gives back the original exactly.

use std::collections::{HashMap, HashSet};

use sbml_io::miriam::{CvTerm, read_annotation, write_annotation};
use sbml_io::{
    Compartment, Document, Event, FunctionDefinition, InitialAssignment, KineticLaw,
    LocalParameter, Meta, Model, ModifierSpeciesReference, Parameter, Reaction, Rule, Species,
    SpeciesReference,
};

/// `orig` with the change from `base` to `mine` applied.
pub fn merge(orig: &Document, base: &Model, mine: &Model) -> Document {
    let mut doc = orig.clone();
    let Some(m) = doc.model.as_mut() else {
        return doc;
    };
    // Symbols whose value is computed (an initial assignment or an
    // assignment rule). Antimony works out a stored value for them from the
    // other values, so editing one value changes theirs; but the stored
    // value is not what the model uses, so it is never copied over.
    let derived: std::collections::HashSet<String> = mine
        .initial_assignments
        .iter()
        .map(|a| a.symbol.clone())
        .chain(mine.rules.iter().filter_map(|r| match r {
            Rule::Assignment { variable, .. } => Some(variable.clone()),
            _ => None,
        }))
        .collect();
    let derived = &derived;
    merge_list(&mut m.function_definitions, &base.function_definitions, &mine.function_definitions, |f, _| f.id.clone(), update_function);
    merge_list(&mut m.unit_definitions, &base.unit_definitions, &mine.unit_definitions, |u, _| u.id.clone(), |o, b, n| {
        set(&mut o.units, &b.units, &n.units);
    });
    merge_list(&mut m.compartments, &base.compartments, &mine.compartments, |c, _| c.id.clone(), |o, b, n| {
        update_compartment(o, b, n, derived.contains(&n.id))
    });
    merge_list(&mut m.species, &base.species, &mine.species, |s, _| s.id.clone(), |o, b, n| {
        update_species(o, b, n, derived.contains(&n.id))
    });
    merge_list(&mut m.parameters, &base.parameters, &mine.parameters, |p, _| p.id.clone(), |o, b, n| {
        update_parameter(o, b, n, derived.contains(&n.id))
    });
    merge_list(&mut m.initial_assignments, &base.initial_assignments, &mine.initial_assignments, |a, _| a.symbol.clone(), update_initial_assignment);
    merge_list(&mut m.rules, &base.rules, &mine.rules, rule_key, update_rule);
    merge_list(&mut m.reactions, &base.reactions, &mine.reactions, |r, _| r.id.clone(), update_reaction);
    merge_list(&mut m.events, &base.events, &mine.events, event_key, update_event);
    merge_links(m, base, mine, &orig.namespaces);
    doc
}

/// An element's database links (biology qualifiers), from its annotation.
fn links_of(meta: &Meta, namespaces: &[(Option<String>, String)]) -> Vec<CvTerm> {
    meta.annotation
        .as_deref()
        .map(|x| read_annotation(x, namespaces).terms)
        .unwrap_or_default()
        .into_iter()
        .filter(|t| !t.model && !t.resources.is_empty())
        .collect()
}

/// Applies edits of database links. Where an element's links in the edited
/// text differ from the text as written, its annotation is rebuilt with the
/// edited links, keeping everything else in it: the model qualifiers, the
/// history, RDF this does not take apart, and other tools' annotations.
fn merge_links(m: &mut Model, base: &Model, mine: &Model, namespaces: &[(Option<String>, String)]) {
    let mut metaids: HashSet<String> = HashSet::new();
    collect_metaids(m, &mut metaids);
    let mut apply = |id: &str, meta: &mut Meta, b: Option<&Meta>, n: Option<&Meta>| {
        let (Some(b), Some(n)) = (b, n) else { return };
        // Text written by Antimony declares its own namespaces.
        let (bt, nt) = (links_of(b, &[]), links_of(n, &[]));
        if bt == nt {
            return;
        }
        rewrite_links(id, meta, nt, namespaces, &mut metaids);
    };
    for c in &mut m.compartments {
        let b = base.compartments.iter().find(|x| x.id == c.id).map(|x| &x.meta);
        let n = mine.compartments.iter().find(|x| x.id == c.id).map(|x| &x.meta);
        apply(&c.id.clone(), &mut c.meta, b, n);
    }
    for s in &mut m.species {
        let b = base.species.iter().find(|x| x.id == s.id).map(|x| &x.meta);
        let n = mine.species.iter().find(|x| x.id == s.id).map(|x| &x.meta);
        apply(&s.id.clone(), &mut s.meta, b, n);
    }
    for p in &mut m.parameters {
        let b = base.parameters.iter().find(|x| x.id == p.id).map(|x| &x.meta);
        let n = mine.parameters.iter().find(|x| x.id == p.id).map(|x| &x.meta);
        apply(&p.id.clone(), &mut p.meta, b, n);
    }
    for r in &mut m.reactions {
        let b = base.reactions.iter().find(|x| x.id == r.id).map(|x| &x.meta);
        let n = mine.reactions.iter().find(|x| x.id == r.id).map(|x| &x.meta);
        apply(&r.id.clone(), &mut r.meta, b, n);
    }
}

fn collect_metaids(m: &Model, out: &mut HashSet<String>) {
    let mut add = |meta: &Meta| {
        if let Some(id) = &meta.metaid {
            out.insert(id.clone());
        }
    };
    add(&m.meta);
    m.compartments.iter().for_each(|x| add(&x.meta));
    m.species.iter().for_each(|x| add(&x.meta));
    m.parameters.iter().for_each(|x| add(&x.meta));
    m.reactions.iter().for_each(|x| add(&x.meta));
    m.function_definitions.iter().for_each(|x| add(&x.meta));
    m.unit_definitions.iter().for_each(|x| add(&x.meta));
    m.events.iter().for_each(|x| add(&x.meta));
}

/// `meta`'s annotation with its biology links replaced by `links`.
fn rewrite_links(
    id: &str,
    meta: &mut Meta,
    links: Vec<CvTerm>,
    namespaces: &[(Option<String>, String)],
    metaids: &mut HashSet<String>,
) {
    let old = meta
        .annotation
        .as_deref()
        .map(|x| read_annotation(x, namespaces))
        .unwrap_or_default();
    let mut terms: Vec<CvTerm> = old.terms.iter().filter(|t| t.model).cloned().collect();
    terms.extend(links);
    // The RDF is about the element by its metaid; give it one if needed.
    let metaid = match &meta.metaid {
        Some(m) => m.clone(),
        None => {
            let mut m = format!("meta_{id}");
            while metaids.contains(&m) {
                m.push('_');
            }
            metaids.insert(m.clone());
            m
        }
    };
    let rdf = write_annotation(&metaid, &terms, &old.history, &old.descriptions);
    let mut parts: Vec<String> = old.other.iter().map(|o| o.xml.clone()).collect();
    if let Some(r) = rdf {
        meta.metaid = Some(metaid);
        parts.push(r);
    }
    meta.annotation = (!parts.is_empty()).then(|| parts.join("\n"));
}

/// Merges one list of elements keyed by `key` (given the element and its
/// position, for elements without ids).
fn merge_list<T: Clone + PartialEq>(
    orig: &mut Vec<T>,
    base: &[T],
    mine: &[T],
    key: impl Fn(&T, usize) -> String,
    update: impl Fn(&mut T, &T, &T),
) {
    let keyed = |v: &[T]| -> Vec<(String, usize)> {
        v.iter().enumerate().map(|(i, x)| (key(x, i), i)).collect()
    };
    let base_k: HashMap<String, usize> = keyed(base).into_iter().collect();
    let mine_k: HashMap<String, usize> = keyed(mine).into_iter().collect();
    // Deleted: in the text as written, gone from the edited text.
    let deleted: Vec<String> = base_k
        .keys()
        .filter(|k| !mine_k.contains_key(*k))
        .cloned()
        .collect();
    let orig_keys: Vec<String> = keyed(orig).into_iter().map(|(k, _)| k).collect();
    let mut keep: Vec<bool> = orig_keys.iter().map(|k| !deleted.contains(k)).collect();
    // Changed: in both, and different.
    for (k, &bi) in &base_k {
        let Some(&mi) = mine_k.get(k) else { continue };
        if base[bi] == mine[mi] {
            continue;
        }
        match orig_keys.iter().position(|o| o == k) {
            Some(oi) => update(&mut orig[oi], &base[bi], &mine[mi]),
            None => orig.push(mine[mi].clone()),
        }
    }
    let mut i = 0;
    orig.retain(|_| {
        let k = keep.get(i).copied().unwrap_or(true);
        i += 1;
        k
    });
    keep.clear();
    // Added: in the edited text only, in its order.
    for (k, mi) in keyed(mine) {
        if base_k.contains_key(&k) {
            continue;
        }
        let present = orig
            .iter()
            .enumerate()
            .any(|(i, x)| key(x, i) == k);
        if !present {
            orig.push(mine[mi].clone());
        }
    }
}

/// Sets `*o` to `m` if the user changed the field (`b` differs from `m`).
fn set<T: Clone + PartialEq>(o: &mut T, b: &T, m: &T) {
    if b != m {
        *o = m.clone();
    }
}

fn update_function(o: &mut FunctionDefinition, b: &FunctionDefinition, m: &FunctionDefinition) {
    set(&mut o.math, &b.math, &m.math);
}

/// Like [`set`], but NaN equals NaN: a value written as NaN is unchanged.
fn set_value(o: &mut Option<f64>, b: &Option<f64>, m: &Option<f64>) {
    let same = match (b, m) {
        (Some(x), Some(y)) => x == y || (x.is_nan() && y.is_nan()),
        (None, None) => true,
        _ => false,
    };
    if !same {
        *o = *m;
    }
}

fn update_compartment(o: &mut Compartment, b: &Compartment, m: &Compartment, derived: bool) {
    if !derived {
        set_value(&mut o.size, &b.size, &m.size);
    }
    set(&mut o.constant, &b.constant, &m.constant);
    set(&mut o.spatial_dimensions, &b.spatial_dimensions, &m.spatial_dimensions);
    set(&mut o.units, &b.units, &m.units);
    set(&mut o.meta.sbo_term, &b.meta.sbo_term, &m.meta.sbo_term);
}

fn update_species(o: &mut Species, b: &Species, m: &Species, derived: bool) {
    set(&mut o.compartment, &b.compartment, &m.compartment);
    // A value is one of amount or concentration: take the pair together.
    let (mut a, mut c) = (o.initial_amount, o.initial_concentration);
    set_value(&mut a, &b.initial_amount, &m.initial_amount);
    set_value(&mut c, &b.initial_concentration, &m.initial_concentration);
    if !derived && (a, c) != (o.initial_amount, o.initial_concentration) {
        o.initial_amount = m.initial_amount;
        o.initial_concentration = m.initial_concentration;
    }
    set(&mut o.boundary_condition, &b.boundary_condition, &m.boundary_condition);
    set(&mut o.constant, &b.constant, &m.constant);
    set(&mut o.has_only_substance_units, &b.has_only_substance_units, &m.has_only_substance_units);
    set(&mut o.substance_units, &b.substance_units, &m.substance_units);
    set(&mut o.meta.sbo_term, &b.meta.sbo_term, &m.meta.sbo_term);
}

fn update_parameter(o: &mut Parameter, b: &Parameter, m: &Parameter, derived: bool) {
    if !derived {
        set_value(&mut o.value, &b.value, &m.value);
    }
    set(&mut o.constant, &b.constant, &m.constant);
    set(&mut o.units, &b.units, &m.units);
    set(&mut o.meta.sbo_term, &b.meta.sbo_term, &m.meta.sbo_term);
}

fn update_initial_assignment(o: &mut InitialAssignment, b: &InitialAssignment, m: &InitialAssignment) {
    set(&mut o.math, &b.math, &m.math);
}

fn rule_key(r: &Rule, i: usize) -> String {
    match r {
        Rule::Assignment { variable, .. } => format!("assignment:{variable}"),
        Rule::Rate { variable, .. } => format!("rate:{variable}"),
        Rule::Algebraic { .. } => format!("algebraic:{i}"),
    }
}

fn update_rule(o: &mut Rule, b: &Rule, m: &Rule) {
    if b.math() == m.math() {
        return;
    }
    let new = m.math().cloned();
    match o {
        Rule::Assignment { math, .. } | Rule::Rate { math, .. } | Rule::Algebraic { math, .. } => {
            *math = new;
        }
    }
}

fn update_reaction(o: &mut Reaction, b: &Reaction, m: &Reaction) {
    set(&mut o.reversible, &b.reversible, &m.reversible);
    set(&mut o.meta.sbo_term, &b.meta.sbo_term, &m.meta.sbo_term);
    if refs_differ(&b.reactants, &m.reactants) {
        o.reactants = merge_refs(&o.reactants, &m.reactants);
    }
    if refs_differ(&b.products, &m.products) {
        o.products = merge_refs(&o.products, &m.products);
    }
    let modifiers = |v: &[ModifierSpeciesReference]| -> Vec<String> {
        v.iter().map(|x| x.species.clone()).collect()
    };
    if modifiers(&b.modifiers) != modifiers(&m.modifiers) {
        o.modifiers = m
            .modifiers
            .iter()
            .map(|mm| {
                o.modifiers
                    .iter()
                    .find(|x| x.species == mm.species)
                    .cloned()
                    .unwrap_or_else(|| mm.clone())
            })
            .collect();
    }
    match (&mut o.kinetic_law, &b.kinetic_law, &m.kinetic_law) {
        (Some(ok), Some(bk), Some(mk)) => update_kinetic_law(ok, bk, mk),
        (Some(ok), None, Some(mk)) => update_kinetic_law(ok, &KineticLaw::default(), mk),
        (None, _, Some(mk)) => o.kinetic_law = Some(mk.clone()),
        (Some(_), Some(_), None) => o.kinetic_law = None,
        _ => {}
    }
}

fn update_kinetic_law(o: &mut KineticLaw, b: &KineticLaw, m: &KineticLaw) {
    set(&mut o.math, &b.math, &m.math);
    merge_list(&mut o.parameters, &b.parameters, &m.parameters, |p, _| p.id.clone(), update_local);
}

fn update_local(o: &mut LocalParameter, b: &LocalParameter, m: &LocalParameter) {
    set_value(&mut o.value, &b.value, &m.value);
}

fn refs_differ(b: &[SpeciesReference], m: &[SpeciesReference]) -> bool {
    let sig = |v: &[SpeciesReference]| -> Vec<(String, Option<f64>)> {
        v.iter().map(|x| (x.species.clone(), x.stoichiometry)).collect()
    };
    sig(b) != sig(m)
}

/// The edited participants, keeping each one's original element (its SBO
/// term, id, annotation) where the species was already there.
fn merge_refs(orig: &[SpeciesReference], mine: &[SpeciesReference]) -> Vec<SpeciesReference> {
    mine.iter()
        .map(|m| match orig.iter().find(|o| o.species == m.species) {
            Some(o) => {
                let mut r = o.clone();
                if r.stoichiometry != m.stoichiometry {
                    r.stoichiometry = m.stoichiometry;
                    r.stoichiometry_math = None;
                }
                r
            }
            None => m.clone(),
        })
        .collect()
}

/// Events are written with their id, or `_E<position>` if they have none.
fn event_key(e: &Event, i: usize) -> String {
    e.id.clone().unwrap_or_else(|| format!("_E{i}"))
}

fn update_event(o: &mut Event, b: &Event, m: &Event) {
    let trig = |e: &Event| e.trigger.as_ref().and_then(|t| t.math.clone());
    if trig(b) != trig(m) {
        let t = o.trigger.get_or_insert_with(Default::default);
        t.math = trig(m);
    }
    let delay = |e: &Event| e.delay.as_ref().and_then(|d| d.math.clone());
    if delay(b) != delay(m) {
        match delay(m) {
            Some(d) => o.delay.get_or_insert_with(Default::default).math = Some(d),
            None => o.delay = None,
        }
    }
    let sets = |e: &Event| -> Vec<(String, Option<sbml_io::Math>)> {
        e.assignments.iter().map(|a| (a.variable.clone(), a.math.clone())).collect()
    };
    if sets(b) != sets(m) {
        o.assignments = m
            .assignments
            .iter()
            .map(|ma| match o.assignments.iter().find(|x| x.variable == ma.variable) {
                Some(x) => {
                    let mut x = x.clone();
                    x.math = ma.math.clone();
                    x
                }
                None => ma.clone(),
            })
            .collect();
    }
}
