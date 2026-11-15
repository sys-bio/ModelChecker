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
    let mut failures = 0;
    for f in &files {
        let Ok(xml) = std::fs::read_to_string(f) else { continue };
        let Ok(read) = sbml_io::read(&xml) else { continue };
        let orig = read.document;
        let Some(model) = &orig.model else { continue };
        models += 1;
        let name = f.file_name().unwrap().to_string_lossy().to_string();
        let text = to_antimony(model).text;

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
}
