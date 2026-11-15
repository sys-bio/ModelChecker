//! A model written out as Antimony, for reading alongside the findings.
//!
//! The text aims to be close to what Tellurium's `getAntimony` gives, and
//! readable first: rate laws keep their local parameters' names (Antimony
//! has no local parameters, so their values are listed in a comment after
//! the reaction). Units, constraints, and annotations are left out. Each
//! element's line is recorded, so a finding can point into the text.

use std::collections::{HashMap, HashSet};

use sbml_io::{Math, Model, Rule};

use crate::fmt::infix;
use crate::sample::fmt_num;

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

/// A copy of the model with ids that are Antimony keywords renamed (with a
/// trailing `_`), and the map back from new ids to the original ones.
fn renamed(model: &Model) -> (Model, HashMap<String, String>) {
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
    let back = map.iter().map(|(k, v)| (v.clone(), k.clone())).collect();
    if map.is_empty() {
        return (model.clone(), back);
    }
    let r = |s: &mut String| {
        if let Some(n) = map.get(s.as_str()) {
            *s = n.clone();
        }
    };
    let rm = |m: &mut Option<Math>| {
        if let Some(x) = m {
            *x = rename_math(x, &map);
        }
    };
    let mut m = model.clone();
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
    (m, back)
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

pub fn to_antimony(model: &Model) -> AntimonyText {
    let (model, back) = renamed(model);
    let mut text = write(&model);
    // Line keys use the model's own ids.
    for (_, id, _) in &mut text.lines {
        if let Some(orig) = back.get(id.as_str()) {
            *id = orig.clone();
        } else if let Some((r, p)) = id.split_once('.') {
            let r = back.get(r).map(String::as_str).unwrap_or(r);
            let p = back.get(p).map(String::as_str).unwrap_or(p);
            *id = format!("{r}.{p}");
        }
    }
    text
}

fn write(model: &Model) -> AntimonyText {
    let mut o = Out {
        text: String::new(),
        line: 0,
        lines: Vec::new(),
    };
    let id = model.id.as_deref().unwrap_or("model");
    if let Some(name) = &model.name {
        o.push(&format!("// {name}"));
    }
    o.push("// Written from SBML by Model Checker; units and annotations are left out.");
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
                            Some(n) if n != 1.0 => format!("{} {name}", fmt_num(n)),
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
            let locals: Vec<String> = r
                .kinetic_law
                .iter()
                .flat_map(|k| &k.parameters)
                .map(|p| {
                    let v = p.value.map(fmt_num).unwrap_or_else(|| "?".into());
                    format!("{} = {v}", p.id)
                })
                .collect();
            if !locals.is_empty() {
                let line = format!("    // local to {}: {}", r.id, locals.join(", "));
                // Each local parameter points at this line.
                for p in r.kinetic_law.iter().flat_map(|k| &k.parameters) {
                    o.lines
                        .push(("local parameter", format!("{}.{}", r.id, p.id), o.line + 1));
                }
                o.push(&line);
            }
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
                (Some(c), _) => fmt_num(c),
                (None, Some(a)) if s.has_only_substance_units == Some(true) => fmt_num(a),
                (None, Some(a)) => format!("{}/{}", fmt_num(a), s.compartment),
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
            o.push_for("compartment", &c.id, &format!("  {} = {};", c.id, fmt_num(c.size.unwrap())));
        }
    }
    if !model.parameters.is_empty() {
        o.blank();
        o.push("  // Parameters:");
        for p in &model.parameters {
            let line = match p.value {
                Some(v) => format!("  {} = {};", p.id, fmt_num(v)),
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
