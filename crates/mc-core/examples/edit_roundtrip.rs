//! For every SBML model in a folder, checks that editing its Antimony text
//! never damages the model:
//!
//! 1. the unedited text merges back to exactly the original document;
//! 2. changing one parameter's value changes that value and nothing else.
//!
//! cargo run --release -p mc-core --features antimony --example edit_roundtrip -- FOLDER

use mc_core::antimony_out::{exact, to_antimony};
use mc_core::from_antimony::apply_edit;

/// The first difference between two models, for reporting.
fn first_difference(a: &sbml_io::Model, b: &sbml_io::Model) -> String {
    macro_rules! cmp {
        ($f:ident) => {
            if a.$f != b.$f {
                let n = a.$f.len().max(b.$f.len());
                for i in 0..n {
                    if a.$f.get(i) != b.$f.get(i) {
                        return format!(
                            "{}[{i}]:\n  original: {:?}\n  merged:   {:?}",
                            stringify!($f),
                            a.$f.get(i),
                            b.$f.get(i)
                        );
                    }
                }
            }
        };
    }
    cmp!(function_definitions);
    cmp!(unit_definitions);
    cmp!(compartments);
    cmp!(species);
    cmp!(parameters);
    cmp!(initial_assignments);
    cmp!(rules);
    cmp!(constraints);
    cmp!(reactions);
    cmp!(events);
    "elsewhere in the model".into()
}

fn main() {
    let dir = std::env::args().nth(1).expect("usage: edit_roundtrip FOLDER");
    let mut files: Vec<_> = std::fs::read_dir(&dir)
        .expect("cannot read folder")
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|e| e == "xml"))
        .collect();
    files.sort();
    let (mut models, mut identity_ok, mut edit_ok, mut edit_tried) = (0, 0, 0, 0);
    let (mut law_ok, mut law_tried) = (0, 0);
    let (mut sbo_ok, mut sbo_tried, mut units_ok, mut units_tried) = (0, 0, 0, 0);
    let (mut links_ok, mut links_tried) = (0, 0);
    let mut failures = 0;
    for f in &files {
        let Ok(xml) = std::fs::read_to_string(f) else { continue };
        let Ok(read) = sbml_io::read(&xml) else { continue };
        let orig = read.document;
        let Some(model) = &orig.model else { continue };
        models += 1;
        let name = f.file_name().unwrap().to_string_lossy().to_string();
        let text = to_antimony(model, &orig.namespaces).text;

        // 1. Unedited.
        match apply_edit(&xml, &text) {
            Ok(e) if sbml_io::write(&e.document) == sbml_io::write(&orig) => identity_ok += 1,
            Ok(e) => {
                failures += 1;
                if failures <= 10 {
                    println!("{name}: unedited text changed the model: {}",
                        first_difference(model, e.document.model.as_ref().unwrap()));
                }
            }
            Err(errs) => {
                failures += 1;
                if failures <= 10 {
                    println!("{name}: unedited text gave errors: {:?}", errs.first().map(|e| &e.message));
                }
            }
        }

        // 4 and 5. One SBO term, and one variable's units, changed.
        // The element written as `id` (keyword ids carry a trailing `_`):
        // its units, if it has a units field, and its annotation block.
        fn find_meta<'m>(
            m: &'m mut sbml_io::Model,
            id: &str,
        ) -> Option<(Option<&'m mut Option<String>>, &'m mut sbml_io::Meta)> {
            let plain = id.trim_end_matches('_');
            let is = |x: &str| x == id || x == plain;
            if let Some(c) = m.compartments.iter_mut().find(|x| is(&x.id)) {
                return Some((Some(&mut c.units), &mut c.meta));
            }
            if let Some(s) = m.species.iter_mut().find(|x| is(&x.id)) {
                return Some((Some(&mut s.substance_units), &mut s.meta));
            }
            if let Some(p) = m.parameters.iter_mut().find(|x| is(&x.id)) {
                return Some((Some(&mut p.units), &mut p.meta));
            }
            if let Some(r) = m.reactions.iter_mut().find(|x| is(&x.id)) {
                return Some((None, &mut r.meta));
            }
            None
        }
        if let Some(line) = text.lines().find(|l| l.contains(".sboTerm = ") && !l.trim_start().starts_with("//")) {
            let (id, rest) = line.trim_start().split_once(".sboTerm = ").unwrap();
            let n: u32 = rest.split(';').next().unwrap().trim().parse().unwrap();
            let new_line = line.replacen(&format!("= {n};"), &format!("= {};", n + 1), 1);
            let mut want = orig.clone();
            if let Some((_, meta)) = find_meta(want.model.as_mut().unwrap(), id) {
                sbo_tried += 1;
                meta.sbo_term = Some(format!("SBO:{:07}", n + 1));
                match apply_edit(&xml, &text.replacen(line, &new_line, 1)) {
                    Ok(e) if sbml_io::write(&e.document) == sbml_io::write(&want) => sbo_ok += 1,
                    Ok(e) => {
                        failures += 1;
                        if failures <= 10 {
                            println!("{name}: editing the SBO term of {id} changed more: {}",
                                first_difference(want.model.as_ref().unwrap(), e.document.model.as_ref().unwrap()));
                        }
                    }
                    Err(errs) => {
                        failures += 1;
                        if failures <= 10 {
                            println!("{name}: SBO edit gave errors: {:?}", errs.first().map(|e| &e.message));
                        }
                    }
                }
            }
        }
        if let Some(line) = text.lines().find(|l| l.contains(" has ") && !l.trim_start().starts_with("//")) {
            let (id, unit) = line.trim_start().split_once(" has ").unwrap();
            let unit = unit.trim_end_matches(';');
            let new_unit = if unit == "dimensionless" { "second" } else { "dimensionless" };
            let new_line = line.replacen(&format!(" has {unit};"), &format!(" has {new_unit};"), 1);
            let mut want = orig.clone();
            if let Some((Some(u), _)) = find_meta(want.model.as_mut().unwrap(), id) {
                units_tried += 1;
                *u = Some(new_unit.to_string());
                match apply_edit(&xml, &text.replacen(line, &new_line, 1)) {
                    Ok(e) if sbml_io::write(&e.document) == sbml_io::write(&want) => units_ok += 1,
                    Ok(e) => {
                        failures += 1;
                        if failures <= 10 {
                            println!("{name}: editing the units of {id} changed more: {}",
                                first_difference(want.model.as_ref().unwrap(), e.document.model.as_ref().unwrap()));
                        }
                    }
                    Err(errs) => {
                        failures += 1;
                        if failures <= 10 {
                            println!("{name}: units edit gave errors: {:?}", errs.first().map(|e| &e.message));
                        }
                    }
                }
            }
        }

        // 6. One database link changed: only that element's annotation (and
        // metaid) may differ, and it must hold the new link, not the old.
        if let Some(line) = text.lines().find(|l| l.contains(" identity \"") && !l.trim_start().starts_with("//")) {
            let id = line.trim_start().split(' ').next().unwrap().to_string();
            let old_uri = line.split('"').nth(1).unwrap().to_string();
            let new_uri = if old_uri.contains("CHEBI:15422") {
                "https://identifiers.org/CHEBI:17234"
            } else {
                "https://identifiers.org/CHEBI:15422"
            };
            let new_line = line.replacen(&old_uri, new_uri, 1);
            match apply_edit(&xml, &text.replacen(line, &new_line, 1)) {
                Ok(e) => {
                    links_tried += 1;
                    let mut got = e.document.clone();
                    let gm = got.model.as_mut().unwrap();
                    let new_meta = find_meta(gm, &id).map(|(_, meta)| meta.clone());
                    let mut copy = orig.clone();
                    let orig_meta = find_meta(copy.model.as_mut().unwrap(), &id).map(|(_, meta)| meta.clone());
                    // Judge by the SBML links themselves: other tools' own
                    // copies (COPASI keeps one) are rightly left alone.
                    let ok_link = new_meta.as_ref().and_then(|m| m.annotation.as_deref()).is_some_and(|a| {
                        let terms = sbml_io::miriam::read_annotation(a, &orig.namespaces).terms;
                        let is: Vec<&String> = terms
                            .iter()
                            .filter(|t| !t.model && t.qualifier == "is")
                            .flat_map(|t| &t.resources)
                            .collect();
                        is.iter().any(|r| *r == new_uri)
                            && (!is.iter().any(|r| **r == old_uri) || line.matches(&old_uri).count() > 1)
                    });
                    // Put the original annotation back: the rest must be untouched.
                    if let (Some((_, meta)), Some(o)) = (find_meta(gm, &id), orig_meta) {
                        *meta = o;
                    }
                    if ok_link && sbml_io::write(&got) == sbml_io::write(&orig) {
                        links_ok += 1;
                    } else {
                        failures += 1;
                        if failures <= 10 {
                            println!("{name}: editing a link of {id} {}", if ok_link {
                                format!("changed more: {}", first_difference(model, got.model.as_ref().unwrap()))
                            } else {
                                "did not change the link".into()
                            });
                        }
                    }
                }
                Err(errs) => {
                    failures += 1;
                    if failures <= 10 {
                        println!("{name}: link edit gave errors: {:?}", errs.first().map(|e| &e.message));
                    }
                }
            }
        }

        // 3. One rate law changed: only that rate law may differ.
        if let Some(r) = model.reactions.iter().find(|r| r.kinetic_law.as_ref().is_some_and(|k| k.math.is_some())) {
            // Keyword ids are written with a trailing `_`.
            let written = text
                .lines()
                .find(|l| {
                    let t = l.trim_start();
                    t.starts_with(&format!("{}: ", r.id)) || t.starts_with(&format!("{}_: ", r.id))
                })
                .map(str::to_string);
            if let Some(line) = written
                && let Some(semi) = line.find("; ")
                && line.ends_with(';')
            {
                law_tried += 1;
                let law = &line[semi + 2..line.len() - 1];
                let new_line = format!("{}; ({law})*2;", &line[..semi]);
                let edited = text.replacen(&line, &new_line, 1);
                match apply_edit(&xml, &edited) {
                    Ok(e) => {
                        let mut got = e.document.clone();
                        let gm = got.model.as_mut().unwrap();
                        let gr = gm.reactions.iter_mut().find(|x| x.id == r.id).unwrap();
                        let changed = gr.kinetic_law.as_ref().unwrap().math != r.kinetic_law.as_ref().unwrap().math;
                        // Put the original rate law back: the rest must be untouched.
                        gr.kinetic_law.as_mut().unwrap().math = r.kinetic_law.as_ref().unwrap().math.clone();
                        if changed && sbml_io::write(&got) == sbml_io::write(&orig) {
                            law_ok += 1;
                        } else {
                            failures += 1;
                            if failures <= 10 {
                                println!("{name}: editing the rate law of {} {}", r.id,
                                    if changed { first_difference(model, got.model.as_ref().unwrap()) } else { "had no effect".into() });
                            }
                        }
                    }
                    Err(errs) => {
                        failures += 1;
                        if failures <= 10 {
                            println!("{name}: edited rate law gave errors: {:?}", errs.first().map(|e| &e.message));
                        }
                    }
                }
            }
        }

        // 2. One parameter's value changed: one whose value is used, that is
        // not computed by an initial assignment or an assignment rule.
        let computed = |id: &str| {
            model.initial_assignments.iter().any(|a| a.symbol == id)
                || model.rules.iter().any(|r| matches!(r, sbml_io::Rule::Assignment { variable, .. } if variable == id))
        };
        let Some((i, p, v)) = model
            .parameters
            .iter()
            .enumerate()
            .filter(|(_, p)| !computed(&p.id))
            .find_map(|(i, p)| p.value.map(|v| (i, p, v)))
        else {
            continue;
        };
        let line = format!("  {} = {};", p.id, exact(v));
        if !text.contains(&format!("\n{line}\n")) {
            continue;
        }
        edit_tried += 1;
        let new_v = if v == 0.0 { 1.0 } else { v * 2.0 };
        let edited = text.replacen(&format!("\n{line}\n"), &format!("\n  {} = {};\n", p.id, exact(new_v)), 1);
        let mut want = orig.clone();
        want.model.as_mut().unwrap().parameters[i].value = Some(new_v);
        match apply_edit(&xml, &edited) {
            Ok(e) if sbml_io::write(&e.document) == sbml_io::write(&want) => edit_ok += 1,
            Ok(e) => {
                failures += 1;
                if failures <= 10 {
                    println!("{name}: editing {} changed more than its value: {}", p.id,
                        first_difference(want.model.as_ref().unwrap(), e.document.model.as_ref().unwrap()));
                }
            }
            Err(errs) => {
                failures += 1;
                if failures <= 10 {
                    println!("{name}: edited text gave errors: {:?}", errs.first().map(|e| &e.message));
                }
            }
        }
    }
    println!("{models} models");
    println!("unedited text gives back the original exactly: {identity_ok} of {models}");
    println!("changing one value changes only that value:   {edit_ok} of {edit_tried}");
    println!("changing one rate law changes only that law:  {law_ok} of {law_tried}");
    println!("changing one SBO term changes only that term: {sbo_ok} of {sbo_tried}");
    println!("changing one unit changes only that unit:     {units_ok} of {units_tried}");
    println!("changing one link changes only that link:     {links_ok} of {links_tried}");
}
