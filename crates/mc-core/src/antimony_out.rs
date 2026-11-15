//! A model written out as Antimony, to read and edit beside the findings.
//!
//! The text is close to what Tellurium's `getAntimony` gives. Two things
//! SBML has and Antimony lacks are handled by a reversible [`prepare`] step:
//! a reaction's local parameter `k1` becomes a global `J1_k1` (as Tellurium
//! writes it), and an id that is an Antimony keyword (a compartment called
//! `compartment`) gets a trailing `_`. [`unprepare`] undoes both, so an
//! edited text can be merged back into the original model without loss
//! (see [`crate::merge`]). Units, SBO terms and annotations are not written;
//! the merge keeps them from the original. Numbers are written exactly.
//! Each element's line is recorded, so a finding can point into the text.

use std::collections::{HashMap, HashSet};

use sbml_io::{Math, Model, Rule};

use crate::fmt::infix;

/// Antimony text, and where each element is in it.
#[derive(Debug, Clone, PartialEq)]
pub struct AntimonyText {
    pub text: String,
    /// (element kind, element id, line from 1). Kinds and ids are those of
    /// [`crate::Element`], so a finding can be looked up.
    pub lines: Vec<(&'static str, String, usize)>,
}

struct Out {
    text: String,
    line: usize,
    lines: Vec<(&'static str, String, usize)>,
}

impl Out {
    fn push(&mut self, s: &str) {
        self.text.push_str(s);
        self.text.push('\n');
        self.line += 1;
    }
    /// A line belonging to an element.
    fn push_for(&mut self, kind: &'static str, id: &str, s: &str) {
        self.lines.push((kind, id.to_string(), self.line + 1));
        self.push(s);
    }
    fn blank(&mut self) {
        if !self.text.ends_with("\n\n") && !self.text.is_empty() {
            self.push("");
        }
    }
}

fn body(m: &Math) -> &Math {
    match m {
        Math::Semantics { math, .. } => body(math),
        _ => m,
    }
}

/// Words Antimony reserves, which an SBML id may still be. `compartment`
/// is a common one: several tools give it as the default compartment's id.
const RESERVED: &[&str] = &[
    "DNA", "after", "at", "avogadro", "compartment", "const", "delete", "end", "event", "ext",
    "formula", "fromTrigger", "function", "gene", "has", "import", "in", "inf", "infinity", "is",
    "model", "module", "nan", "notanumber", "operator", "persistent", "priority", "reaction",
    "species", "substanceOnly", "t0", "time", "unit", "var",
];

/// A number written so that it reads back as exactly the same value.
pub fn exact(x: f64) -> String {
    if x == 0.0 {
        "0".into()
    } else if x.fract() == 0.0 && x.abs() < 1e15 {
        format!("{}", x as i64)
    } else if (1e-4..1e15).contains(&x.abs()) {
        // Rust prints the shortest text that reads back as the same f64.
        format!("{x}")
    } else {
        format!("{x:e}")
    }
}

/// A local parameter written as a global: `local` in reaction `reaction`
/// is the global `global` in the text.
#[derive(Debug, Clone, PartialEq)]
pub struct Local {
    pub reaction: String,
    pub local: String,
    pub global: String,
}

/// How a model was prepared for Antimony, so it can be undone.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Mapping {
    /// Keyword renames: new id to original id.
    pub back: HashMap<String, String>,
    pub locals: Vec<Local>,
}

impl Mapping {
    /// The element a written (kind, id) stands for in the original model.
    pub fn key(&self, kind: &'static str, id: &str) -> (&'static str, String) {
        let id = self.back.get(id).map(String::as_str).unwrap_or(id);
        if kind == "parameter"
            && let Some(l) = self.locals.iter().find(|l| l.global == id)
        {
            return ("local parameter", format!("{}.{}", l.reaction, l.local));
        }
        (kind, id.to_string())
    }
}

/// The model ready to write as Antimony: local parameters made global and
/// keyword ids renamed. The mapping undoes both ([`unprepare`]).
pub fn prepare(model: &Model) -> (Model, Mapping) {
    let (mut m, locals) = globalize(model);
    let map = keyword_map(&m);
    rename_model(&mut m, &map);
    let back = map.into_iter().map(|(k, v)| (v, k)).collect();
    (m, Mapping { back, locals })
}

/// Undoes [`prepare`] on a model read back from (possibly edited) Antimony.
pub fn unprepare(m: &mut Model, map: &Mapping) {
    rename_model(m, &map.back);
    localize(m, &map.locals);
}

/// Every id in the model, for avoiding clashes.
fn all_ids(model: &Model) -> HashSet<String> {
    let mut ids: HashSet<String> = HashSet::new();
    ids.extend(model.function_definitions.iter().map(|f| f.id.clone()));
    ids.extend(model.compartments.iter().map(|c| c.id.clone()));
    ids.extend(model.species.iter().map(|s| s.id.clone()));
    ids.extend(model.parameters.iter().map(|p| p.id.clone()));
    ids.extend(model.reactions.iter().map(|r| r.id.clone()));
    ids.extend(model.events.iter().filter_map(|e| e.id.clone()));
    for r in &model.reactions {
        ids.extend(r.kinetic_law.iter().flat_map(|k| &k.parameters).map(|p| p.id.clone()));
    }
    ids
}

/// Local parameters as globals named `reaction_local`.
fn globalize(model: &Model) -> (Model, Vec<Local>) {
    let mut m = model.clone();
    let mut used = all_ids(model);
    let mut locals = Vec::new();
    for r in &mut m.reactions {
        let Some(k) = &mut r.kinetic_law else { continue };
        let mut map: HashMap<String, String> = HashMap::new();
        for p in std::mem::take(&mut k.parameters) {
            let mut global = format!("{}_{}", r.id, p.id);
            while used.contains(&global) {
                global.push('_');
            }
            used.insert(global.clone());
            m.parameters.push(sbml_io::Parameter {
                id: global.clone(),
                name: p.name.clone(),
                value: p.value,
                units: p.units.clone(),
                constant: Some(true),
                meta: Default::default(),
            });
            map.insert(p.id.clone(), global.clone());
            locals.push(Local {
                reaction: r.id.clone(),
                local: p.id,
                global,
            });
        }
        if let Some(math) = &mut k.math {
            *math = rename_math(math, &map);
        }
    }
    (m, locals)
}

/// Globals made from local parameters turned back into local parameters.
fn localize(m: &mut Model, locals: &[Local]) {
    for l in locals {
        let Some(i) = m.parameters.iter().position(|p| p.id == l.global) else {
            continue;
        };
        let p = m.parameters.remove(i);
        let Some(r) = m.reactions.iter_mut().find(|r| r.id == l.reaction) else {
            continue;
        };
        let k = r.kinetic_law.get_or_insert_with(Default::default);
        if let Some(math) = &mut k.math {
            let map = HashMap::from([(l.global.clone(), l.local.clone())]);
            *math = rename_math(math, &map);
        }
        k.parameters.push(sbml_io::LocalParameter {
            id: l.local.clone(),
            name: p.name,
            value: p.value,
            units: p.units,
            meta: Default::default(),
        });
    }
}

/// Renames for ids that are Antimony keywords: original to new.
fn keyword_map(model: &Model) -> HashMap<String, String> {
    let mut used: HashSet<String> = HashSet::new();
    let mut ids: Vec<&str> = Vec::new();
    ids.extend(model.function_definitions.iter().map(|f| f.id.as_str()));
    ids.extend(model.compartments.iter().map(|c| c.id.as_str()));
    ids.extend(model.species.iter().map(|s| s.id.as_str()));
    ids.extend(model.parameters.iter().map(|p| p.id.as_str()));
    ids.extend(model.reactions.iter().map(|r| r.id.as_str()));
    for r in &model.reactions {
        ids.extend(r.kinetic_law.iter().flat_map(|k| &k.parameters).map(|p| p.id.as_str()));
    }
    ids.extend(model.events.iter().filter_map(|e| e.id.as_deref()));
    // Function arguments too: they are local, but are written out.
    for f in &model.function_definitions {
        if let Some(Math::Lambda { args, .. }) = f.math.as_ref().map(body) {
            ids.extend(args.iter().map(String::as_str));
        }
    }
    used.extend(ids.iter().map(|s| s.to_string()));
    let mut map: HashMap<String, String> = HashMap::new();
    for id in ids {
        if RESERVED.contains(&id) && !map.contains_key(id) {
            let mut new = format!("{id}_");
            while used.contains(&new) {
                new.push('_');
            }
            used.insert(new.clone());
            map.insert(id.to_string(), new);
        }
    }
    map
}

/// Renames ids everywhere in the model, as declared and as used.
fn rename_model(m: &mut Model, map: &HashMap<String, String>) {
    if map.is_empty() {
        return;
    }
    let r = |s: &mut String| {
        if let Some(n) = map.get(s.as_str()) {
            *s = n.clone();
        }
    };
    let rm = |m: &mut Option<Math>| {
        if let Some(x) = m {
            *x = rename_math(x, map);
        }
    };
    for f in &mut m.function_definitions {
        r(&mut f.id);
        rm(&mut f.math);
    }
    for c in &mut m.compartments {
        r(&mut c.id);
    }
    for s in &mut m.species {
        r(&mut s.id);
        r(&mut s.compartment);
    }
    for p in &mut m.parameters {
        r(&mut p.id);
    }
    for ia in &mut m.initial_assignments {
        r(&mut ia.symbol);
        rm(&mut ia.math);
    }
    for rule in &mut m.rules {
        match rule {
            Rule::Assignment { variable, math, .. } | Rule::Rate { variable, math, .. } => {
                r(variable);
                rm(math);
            }
            Rule::Algebraic { math, .. } => rm(math),
        }
    }
    for re in &mut m.reactions {
        r(&mut re.id);
        for sr in re.reactants.iter_mut().chain(re.products.iter_mut()) {
            r(&mut sr.species);
        }
        for md in &mut re.modifiers {
            r(&mut md.species);
        }
        if let Some(k) = &mut re.kinetic_law {
            for p in &mut k.parameters {
                r(&mut p.id);
            }
            rm(&mut k.math);
        }
    }
    for e in &mut m.events {
        if let Some(id) = &mut e.id {
            r(id);
        }
        if let Some(t) = &mut e.trigger {
            rm(&mut t.math);
        }
        if let Some(d) = &mut e.delay {
            rm(&mut d.math);
        }
        for a in &mut e.assignments {
            r(&mut a.variable);
            rm(&mut a.math);
        }
    }
}

fn rename_math(m: &Math, map: &HashMap<String, String>) -> Math {
    let rec = |x: &Math| rename_math(x, map);
    let name = |s: &String| map.get(s).cloned().unwrap_or_else(|| s.clone());
    match m {
        Math::Ident(id) => Math::Ident(name(id)),
        Math::Apply(op, a) => Math::Apply(*op, a.iter().map(rec).collect()),
        Math::Call(f, a) => Math::Call(name(f), a.iter().map(rec).collect()),
        Math::CallSymbol(s, a) => Math::CallSymbol(s.clone(), a.iter().map(rec).collect()),
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
        Math::Lambda { args, body } => Math::Lambda {
            args: args.iter().map(name).collect(),
            body: Box::new(rec(body)),
        },
        Math::Semantics { math, annotations } => Math::Semantics {
            math: Box::new(rec(math)),
            annotations: annotations.clone(),
        },
        _ => m.clone(),
    }
}

/// The model as Antimony, with each element's line keyed by the model's own
/// kinds and ids.
pub fn to_antimony(model: &Model) -> AntimonyText {
    let (prepared, map) = prepare(model);
    let mut text = write(&prepared);
    for entry in &mut text.lines {
        let (kind, id) = map.key(entry.0, &entry.1);
        entry.0 = kind;
        entry.1 = id;
    }
    text
}

/// A prepared model ([`prepare`]) as Antimony; line keys are the prepared
/// model's ids.
pub fn write(model: &Model) -> AntimonyText {
    let mut o = Out {
        text: String::new(),
        line: 0,
        lines: Vec::new(),
    };
    let id = model.id.as_deref().unwrap_or("model");
    if let Some(name) = &model.name {
        o.push(&format!("// {name}"));
    }
    o.push("// Written from SBML by Model Checker. Units, SBO terms and annotations are");
    o.push("// not shown, but are kept when you edit. Local parameters are written as");
    o.push("// reaction_name (J1_k1).");
    o.blank();

    // Function definitions come before the model in Antimony.
    for f in &model.function_definitions {
        let Some(Math::Lambda { args, body: b }) = f.math.as_ref().map(body) else {
            continue;
        };
        o.push_for("function", &f.id, &format!("function {}({})", f.id, args.join(", ")));
        o.push(&format!("  {}", infix(b)));
        o.push("end");
        o.blank();
    }

    o.push(&format!("model *{id}()"));

    // Compartments and species.
    if !model.compartments.is_empty() || !model.species.is_empty() {
        o.blank();
        o.push("  // Compartments and species:");
        for c in &model.compartments {
            let mut s = String::from("  ");
            if c.constant == Some(false) {
                s.push_str("var ");
            }
            s.push_str(&format!("compartment {};", c.id));
            o.push_for("compartment", &c.id, &s);
        }
        for sp in &model.species {
            let mut s = String::from("  ");
            if sp.constant == Some(true) {
                s.push_str("const ");
            }
            if sp.has_only_substance_units == Some(true) {
                s.push_str("substanceOnly ");
            }
            s.push_str("species ");
            if sp.boundary_condition == Some(true) {
                s.push('$');
            }
            s.push_str(&format!("{} in {};", sp.id, sp.compartment));
            o.push_for("species", &sp.id, &s);
        }
    }

    // Reactions.
    if !model.reactions.is_empty() {
        o.blank();
        o.push("  // Reactions:");
        for r in &model.reactions {
            let side = |refs: &[sbml_io::SpeciesReference]| {
                refs.iter()
                    .map(|sr| {
                        let boundary = model
                            .species
                            .iter()
                            .find(|s| s.id == sr.species)
                            .is_some_and(|s| s.boundary_condition == Some(true));
                        let name = if boundary {
                            format!("${}", sr.species)
                        } else {
                            sr.species.clone()
                        };
                        match sr.stoichiometry {
                            Some(n) if n != 1.0 => format!("{} {name}", exact(n)),
                            _ => name,
                        }
                    })
                    .collect::<Vec<_>>()
                    .join(" + ")
            };
            let arrow = if r.reversible == Some(false) { "=>" } else { "->" };
            let lhs = side(&r.reactants);
            let rhs = side(&r.products);
            let law = r
                .kinetic_law
                .as_ref()
                .and_then(|k| k.math.as_ref())
                .map(infix)
                .unwrap_or_default();
            let mut s = format!("  {}: ", r.id);
            if !lhs.is_empty() {
                s.push_str(&lhs);
                s.push(' ');
            }
            // `S1 => S2; k*S1`, and with no products `S4 => ; k*S4`.
            s.push_str(arrow);
            s.push(' ');
            s.push_str(&rhs);
            s.push_str(&format!("; {law};"));
            o.push_for("reaction", &r.id, &s);
        }
    }

    // Rules.
    let assignments: Vec<&Rule> = model
        .rules
        .iter()
        .filter(|r| matches!(r, Rule::Assignment { .. }))
        .collect();
    if !assignments.is_empty() {
        o.blank();
        o.push("  // Assignment rules:");
        for r in assignments {
            if let (Some(v), Some(m)) = (r.variable(), r.math()) {
                o.push_for("rule", v, &format!("  {v} := {};", infix(m)));
            }
        }
    }
    let rates: Vec<&Rule> = model
        .rules
        .iter()
        .filter(|r| matches!(r, Rule::Rate { .. }))
        .collect();
    if !rates.is_empty() {
        o.blank();
        o.push("  // Rate rules:");
        for r in rates {
            if let (Some(v), Some(m)) = (r.variable(), r.math()) {
                o.push_for("rule", v, &format!("  {v}' = {};", infix(m)));
            }
        }
    }
    let algebraic: Vec<(usize, &Rule)> = model
        .rules
        .iter()
        .enumerate()
        .filter(|(_, r)| matches!(r, Rule::Algebraic { .. }))
        .collect();
    if !algebraic.is_empty() {
        o.blank();
        o.push("  // Algebraic rules:");
        for (i, r) in algebraic {
            if let Some(m) = r.math() {
                o.push_for("algebraic rule", &format!("#{}", i + 1), &format!("  0 = {};", infix(m)));
            }
        }
    }

    // Events.
    if !model.events.is_empty() {
        o.blank();
        o.push("  // Events:");
        for (i, e) in model.events.iter().enumerate() {
            let id = e.id.clone().unwrap_or_else(|| format!("#{}", i + 1));
            let label = e.id.clone().unwrap_or_else(|| format!("_E{i}"));
            let trigger = e
                .trigger
                .as_ref()
                .and_then(|t| t.math.as_ref())
                .map(infix)
                // Antimony needs a trigger; an SBML event may lack one.
                .unwrap_or_else(|| "false".to_string());
            let when = match e.delay.as_ref().and_then(|d| d.math.as_ref()) {
                Some(d) => format!("{} after {trigger}", infix(d)),
                None => trigger,
            };
            let sets: Vec<String> = e
                .assignments
                .iter()
                .map(|a| {
                    let v = a.math.as_ref().map(infix).unwrap_or_default();
                    format!("{} = {v}", a.variable)
                })
                .collect();
            o.push_for(
                "event",
                &id,
                &format!("  {label}: at {when}: {};", sets.join(", ")),
            );
        }
    }

    // Values.
    let species_init: Vec<String> = model
        .species
        .iter()
        .filter_map(|s| {
            let v = match (s.initial_concentration, s.initial_amount) {
                (Some(c), _) => exact(c),
                (None, Some(a)) if s.has_only_substance_units == Some(true) => exact(a),
                (None, Some(a)) => format!("{}/{}", exact(a), s.compartment),
                (None, None) => return None,
            };
            Some(format!("  {} = {v};", s.id))
        })
        .collect();
    if !species_init.is_empty() {
        o.blank();
        o.push("  // Species initializations:");
        for (s, line) in model
            .species
            .iter()
            .filter(|s| s.initial_concentration.is_some() || s.initial_amount.is_some())
            .zip(species_init)
        {
            o.push_for("species", &s.id, &line);
        }
    }
    let sized: Vec<&sbml_io::Compartment> =
        model.compartments.iter().filter(|c| c.size.is_some()).collect();
    if !sized.is_empty() {
        o.blank();
        o.push("  // Compartment initializations:");
        for c in sized {
            o.push_for("compartment", &c.id, &format!("  {} = {};", c.id, exact(c.size.unwrap())));
        }
    }
    if !model.parameters.is_empty() {
        o.blank();
        o.push("  // Parameters:");
        for p in &model.parameters {
            let line = match p.value {
                Some(v) => format!("  {} = {};", p.id, exact(v)),
                None => format!("  {};", p.id),
            };
            o.push_for("parameter", &p.id, &line);
        }
    }
    if !model.initial_assignments.is_empty() {
        o.blank();
        o.push("  // Initial assignments:");
        for ia in &model.initial_assignments {
            if let Some(m) = &ia.math {
                o.push_for(
                    "initial assignment",
                    &ia.symbol,
                    &format!("  {} = {};", ia.symbol, infix(m)),
                );
            }
        }
    }
    o.blank();
    o.push("end");
    AntimonyText {
        text: o.text,
        lines: o.lines,
    }
}

impl AntimonyText {
    /// The line an element is on: its first, if it has several (a species'
    /// declaration before its initial value).
    pub fn line_of(&self, kind: &str, id: &str) -> Option<usize> {
        self.lines
            .iter()
            .find(|(k, i, _)| *k == kind && i == id)
            .map(|(_, _, l)| *l)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_and_maps() {
        let xml = r#"<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
  <model id="m">
    <listOfCompartments><compartment id="c" size="1" constant="true"/></listOfCompartments>
    <listOfSpecies>
      <species id="S" compartment="c" initialConcentration="1" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
      <species id="P" compartment="c" initialConcentration="0" hasOnlySubstanceUnits="false" boundaryCondition="true" constant="false"/>
    </listOfSpecies>
    <listOfParameters><parameter id="k" value="0.5" constant="true"/></listOfParameters>
    <listOfReactions>
      <reaction id="J1" reversible="false">
        <listOfReactants><speciesReference species="S" stoichiometry="2" constant="true"/></listOfReactants>
        <listOfProducts><speciesReference species="P" stoichiometry="1" constant="true"/></listOfProducts>
        <kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><apply><times/><ci>k</ci><ci>S</ci></apply></math></kineticLaw>
      </reaction>
    </listOfReactions>
  </model>
</sbml>"#;
        let doc = sbml_io::read(xml).unwrap().document;
        let a = to_antimony(doc.model.as_ref().unwrap());
        assert!(a.text.contains("  J1: 2 S => $P; k*S;"), "{}", a.text);
        assert!(a.text.contains("  species $P in c;"), "{}", a.text);
        assert!(a.text.contains("  k = 0.5;"), "{}", a.text);
        let line = a.line_of("reaction", "J1").unwrap();
        assert_eq!(a.text.lines().nth(line - 1).unwrap(), "  J1: 2 S => $P; k*S;");
    }
}
