//! Checking an edit of a model's Antimony text (the web page's panel).
//!
//! The original SBML stays the master copy: the edit is merged into it
//! ([`crate::merge`]), so whatever the text does not show (units, SBO terms,
//! annotations, notes, names) is kept, and an unedited text gives back the
//! original exactly. Antimony's own errors and warnings on the edited text
//! become findings (AN001, AN002) at their lines, and every other finding is
//! placed on the line where its element is written.

use std::collections::HashSet;

use antimony::Severity as AntSeverity;
use sbml_io::Document;

use crate::antimony_out::{self, AntimonyText, Mapping};
use crate::codes::{AN_ERROR, AN_WARNING, SB_READ};
use crate::{Element, Finding, Report, Severity};

/// The original model with an edit of its Antimony text applied, and where
/// each element is in the edited text.
pub struct Edited {
    pub document: Document,
    pub lines: Vec<(&'static str, String, usize)>,
    /// Antimony's warnings on the edited text (AN002).
    pub warnings: Vec<Finding>,
}

fn finding(code: &'static str, severity: Severity, message: String, line: Option<usize>) -> Finding {
    Finding {
        code,
        severity,
        element: Element::Model,
        message,
        witness: None,
        reaction: None,
        title: None,
        line,
    }
}

/// Applies `edited`, an edit of the Antimony text of the model in
/// `original_xml`, to that model. Fails with findings if the original
/// cannot be read or the edited text has errors.
pub fn apply_edit(original_xml: &str, edited: &str) -> Result<Edited, Vec<Finding>> {
    let orig = match sbml_io::read(original_xml) {
        Ok(r) => r.document,
        Err(diags) => {
            return Err(diags
                .into_iter()
                .map(|d| finding(SB_READ, Severity::Error, d.message, None))
                .collect());
        }
    };
    let Some(orig_model) = &orig.model else {
        return Err(vec![finding(SB_READ, Severity::Error, "The original document has no model.".into(), None)]);
    };
    let (prepared, map) = antimony_out::prepare(orig_model);
    let base_text = antimony_out::write(&prepared).text;

    let mine = antimony::parse(edited);
    let base = antimony::parse(&base_text);
    // Warnings the unedited text also gives come from the conversion, not
    // from the user's edit.
    let base_warnings: HashSet<String> =
        base.diagnostics.iter().map(|d| d.message.clone()).collect();
    let mut errors = Vec::new();
    let mut warnings = Vec::new();
    for d in &mine.diagnostics {
        let line = (!d.span.is_unknown()).then_some(d.span.start.line as usize);
        match d.severity {
            AntSeverity::Error => errors.push(finding(AN_ERROR, Severity::Error, d.message.clone(), line)),
            AntSeverity::Warning if !base_warnings.contains(&d.message) => {
                warnings.push(finding(AN_WARNING, Severity::Warning, d.message.clone(), line))
            }
            _ => {}
        }
    }
    if !errors.is_empty() {
        return Err(errors);
    }
    if base.has_errors() {
        return Err(vec![finding(
            AN_ERROR,
            Severity::Error,
            "The model's own Antimony text could not be read back, so edits cannot be applied safely."
                .into(),
            None,
        )]);
    }
    let mine_module = antimony::flatten::flatten(&mine.document);
    let base_module = antimony::flatten::flatten(&base.document);
    let to_model = |m: &antimony::model::Module| -> sbml_io::Model {
        let mut model = antimony::sbml::to_sbml(m).model.unwrap_or_default();
        antimony_out::unprepare(&mut model, &map);
        model
    };
    let document = crate::merge::merge(&orig, &to_model(&base_module), &to_model(&mine_module));
    Ok(Edited {
        document,
        lines: lines(&mine_module, &map),
        warnings,
    })
}

/// Checks the original model with an edit of its Antimony text applied.
pub fn check_edit(original_xml: &str, edited: &str) -> Report {
    let text = AntimonyText {
        text: edited.to_string(),
        lines: Vec::new(),
    };
    match apply_edit(original_xml, edited) {
        Ok(e) => {
            let mut report = crate::check_document(&e.document);
            report.findings.extend(e.warnings);
            report.antimony = Some(AntimonyText {
                lines: e.lines,
                ..text
            });
            crate::sort(&mut report.findings);
            report
        }
        Err(findings) => {
            let mut report = Report {
                level: None,
                version: None,
                model_id: None,
                model_name: None,
                findings,
                antimony: Some(text),
            };
            crate::sort(&mut report.findings);
            report
        }
    }
}

/// The original model with the edit applied, as SBML; `None` if the edited
/// text has errors.
pub fn edited_sbml(original_xml: &str, edited: &str) -> Option<String> {
    apply_edit(original_xml, edited)
        .ok()
        .map(|e| sbml_io::write(&e.document))
}

/// Where each element is written in the text, keyed by the original
/// model's kinds and ids.
fn lines(m: &antimony::model::Module, map: &Mapping) -> Vec<(&'static str, String, usize)> {
    let mut out = Vec::new();
    let line = |s: antimony::Span| (!s.is_unknown()).then_some(s.start.line as usize);
    let mut push = |kind: &'static str, id: &str, l: usize| {
        let (kind, id) = map.key(kind, id);
        out.push((kind, id, l));
    };
    for s in &m.symbols {
        let Some(l) = line(s.decl_span).or(line(s.span)) else { continue };
        // A symbol becomes one of these in SBML; findings name it by that
        // kind, so every kind points at the same line.
        for kind in ["species", "parameter", "compartment", "rule", "initial assignment"] {
            push(kind, &s.name, l);
        }
    }
    for r in &m.reactions {
        if let Some(l) = line(r.span) {
            push("reaction", &r.name, l);
        }
    }
    for (i, e) in m.events.iter().enumerate() {
        if let Some(l) = line(e.span) {
            let id = e.name.clone().unwrap_or_else(|| format!("#{}", i + 1));
            push("event", &id, l);
        }
    }
    for f in &m.functions {
        if let Some(l) = line(f.span) {
            push("function", &f.name, l);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const XML: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
  <model id="m" substanceUnits="mole" timeUnits="second" extentUnits="mole">
    <listOfUnitDefinitions>
      <unitDefinition id="per_s"><listOfUnits><unit kind="second" exponent="-1" scale="0" multiplier="1"/></listOfUnits></unitDefinition>
    </listOfUnitDefinitions>
    <listOfCompartments><compartment id="compartment" size="1" constant="true"/></listOfCompartments>
    <listOfSpecies>
      <species id="S" name="Substrate" compartment="compartment" initialConcentration="1" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false" sboTerm="SBO:0000247"/>
      <species id="P" compartment="compartment" initialConcentration="0" hasOnlySubstanceUnits="false" boundaryCondition="false" constant="false"/>
    </listOfSpecies>
    <listOfParameters><parameter id="k" value="0.123456789" units="per_s" constant="true"/></listOfParameters>
    <listOfReactions>
      <reaction id="J1" reversible="false" sboTerm="SBO:0000176">
        <listOfReactants><speciesReference species="S" stoichiometry="1" constant="true"/></listOfReactants>
        <listOfProducts><speciesReference species="P" stoichiometry="1" constant="true"/></listOfProducts>
        <kineticLaw sboTerm="SBO:0000049"><math xmlns="http://www.w3.org/1998/Math/MathML"><apply><times/><ci>k</ci><ci>kl</ci><ci>S</ci></apply></math>
          <listOfLocalParameters><localParameter id="kl" value="2"/></listOfLocalParameters>
        </kineticLaw>
      </reaction>
    </listOfReactions>
  </model>
</sbml>"#;

    fn text() -> String {
        let doc = sbml_io::read(XML).unwrap().document;
        antimony_out::to_antimony(doc.model.as_ref().unwrap()).text
    }

    #[test]
    fn unedited_text_gives_back_the_original() {
        let orig = sbml_io::read(XML).unwrap().document;
        let e = apply_edit(XML, &text()).ok().unwrap();
        assert_eq!(e.document, orig);
    }

    #[test]
    fn an_edit_changes_only_what_was_edited() {
        let edited = text().replace("k = 0.123456789;", "k = 0.5;");
        assert_ne!(edited, text());
        let e = apply_edit(XML, &edited).ok().unwrap();
        let mut want = sbml_io::read(XML).unwrap().document;
        want.model.as_mut().unwrap().parameters[0].value = Some(0.5);
        assert_eq!(e.document, want);
        // Units, SBO terms and names are still there.
        let m = e.document.model.as_ref().unwrap();
        assert_eq!(m.parameters[0].units.as_deref(), Some("per_s"));
        assert_eq!(m.species[0].name.as_deref(), Some("Substrate"));
        assert_eq!(m.reactions[0].kinetic_law.as_ref().unwrap().meta.sbo_term.as_deref(), Some("SBO:0000049"));
    }

    #[test]
    fn units_and_sbo_terms_are_shown_and_editable() {
        let t = text();
        for want in [
            "unit per_s = second^-1;",
            "k has per_s;",
            "S.sboTerm = 247;",
            "J1.sboTerm = 176;",
            "// J1 rate law: SBO:0000049",
            "// Model units (not editable here): substance mole, time second, extent mole",
        ] {
            assert!(t.contains(want), "missing {want:?} in:\n{t}");
        }
        // Change one SBO term: only it changes.
        let e = apply_edit(XML, &t.replace("S.sboTerm = 247;", "S.sboTerm = 245;")).ok().unwrap();
        let mut want = sbml_io::read(XML).unwrap().document;
        want.model.as_mut().unwrap().species[0].meta.sbo_term = Some("SBO:0000245".into());
        assert_eq!(sbml_io::write(&e.document), sbml_io::write(&want));
        // Change a parameter's units: only they change.
        let e = apply_edit(XML, &t.replace("k has per_s;", "k has second;")).ok().unwrap();
        let mut want = sbml_io::read(XML).unwrap().document;
        want.model.as_mut().unwrap().parameters[0].units = Some("second".into());
        assert_eq!(sbml_io::write(&e.document), sbml_io::write(&want));
    }

    #[test]
    fn local_parameters_survive_editing() {
        let t = text();
        assert!(t.contains("J1_kl = 2;"), "{t}");
        let edited = t.replace("J1_kl = 2;", "J1_kl = 3;");
        let e = apply_edit(XML, &edited).ok().unwrap();
        let k = e.document.model.as_ref().unwrap().reactions[0].kinetic_law.clone().unwrap();
        assert_eq!(k.parameters.len(), 1);
        assert_eq!((k.parameters[0].id.as_str(), k.parameters[0].value), ("kl", Some(3.0)));
        // The rate law still names the local parameter as `kl`.
        assert!(crate::fmt::infix(k.math.as_ref().unwrap()).contains("kl"));
        assert!(!e.document.model.as_ref().unwrap().parameters.iter().any(|p| p.id == "J1_kl"));
    }

    #[test]
    fn keyword_ids_survive_editing() {
        let t = text();
        assert!(t.contains("compartment compartment_;"), "{t}");
        let edited = t.replace("compartment_ = 1;", "compartment_ = 2;");
        let e = apply_edit(XML, &edited).ok().unwrap();
        let c = &e.document.model.as_ref().unwrap().compartments[0];
        assert_eq!((c.id.as_str(), c.size), ("compartment", Some(2.0)));
    }

    #[test]
    fn adding_and_deleting_reactions() {
        let t = text();
        let added = t.replace("  // Species", "  J2: P => ; k*P;\n\n  // Species");
        let e = apply_edit(XML, &added).ok().unwrap();
        let ids: Vec<&str> = e.document.model.as_ref().unwrap().reactions.iter().map(|r| r.id.as_str()).collect();
        assert_eq!(ids, ["J1", "J2"]);
        let line = t.lines().find(|l| l.trim_start().starts_with("J1:")).unwrap();
        let deleted = t.replace(&format!("{line}\n"), "");
        let e = apply_edit(XML, &deleted).ok().unwrap();
        assert!(e.document.model.as_ref().unwrap().reactions.is_empty());
    }

    #[test]
    fn findings_on_the_users_lines_and_errors() {
        let edited = text().replace("k*J1_kl*S", "k*J1_kl");
        let r = check_edit(XML, &edited);
        let line = edited.lines().position(|l| l.trim_start().starts_with("J1:")).unwrap() + 1;
        let f = r.findings.iter().find(|f| f.code == crate::codes::RL_ZERO_ORDER).expect("zero order");
        let a = r.antimony.as_ref().unwrap();
        assert_eq!(a.line_of(f.element.kind(), &f.element.id()), Some(line));
        let r = check_edit(XML, &text().replace("k*J1_kl*S;", "k*(J1_kl*S;"));
        assert!(r.findings.iter().any(|f| f.code == AN_ERROR && f.line == Some(line)));
    }
}
