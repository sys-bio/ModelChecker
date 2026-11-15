//! Static checks for kinetic SBML models.
//!
//! ```
//! let xml = r#"<?xml version="1.0" encoding="UTF-8"?>
//! <sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
//!   <model id="m">
//!     <listOfCompartments><compartment id="c" size="1" constant="true"/></listOfCompartments>
//!     <listOfSpecies>
//!       <species id="S" compartment="c" initialConcentration="1" hasOnlySubstanceUnits="false"
//!                boundaryCondition="false" constant="false"/>
//!     </listOfSpecies>
//!     <listOfParameters><parameter id="k" value="1" constant="true"/></listOfParameters>
//!     <listOfReactions>
//!       <reaction id="J" reversible="false">
//!         <listOfReactants><speciesReference species="S" stoichiometry="1" constant="true"/></listOfReactants>
//!         <kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><ci>k</ci></math></kineticLaw>
//!       </reaction>
//!     </listOfReactions>
//!   </model>
//! </sbml>"#;
//! let report = mc_core::check(xml);
//! assert!(report.findings.iter().any(|f| f.code == mc_core::codes::RL_ZERO_ORDER));
//! ```

pub mod antimony_out;
pub mod checks;
pub mod codes;
pub mod eval;
pub mod finding;
pub mod fmt;
#[cfg(feature = "antimony")]
pub mod from_antimony;
pub mod json;
pub mod merge;
pub mod moiety;
pub mod sample;
pub mod sbo;
pub mod sbo_data;
pub mod registry;
pub mod registry_data;
pub mod units;
pub mod value;
pub mod view;

pub use finding::{Element, Finding, Severity};

/// The result of checking one document.
#[derive(Debug, Clone, PartialEq)]
pub struct Report {
    pub level: Option<u32>,
    pub version: Option<u32>,
    pub model_id: Option<String>,
    pub model_name: Option<String>,
    pub findings: Vec<Finding>,
    /// The model as Antimony, with each element's line, for reading
    /// alongside the findings.
    pub antimony: Option<antimony_out::AntimonyText>,
}

impl Report {
    pub fn count(&self, severity: Severity) -> usize {
        self.findings.iter().filter(|f| f.severity == severity).count()
    }
}

/// Reads an SBML document and runs every check on it. Never fails: a file
/// that cannot be read gives a report with one SB001 finding.
pub fn check(xml: &str) -> Report {
    let mut report = Report {
        level: None,
        version: None,
        model_id: None,
        model_name: None,
        findings: Vec::new(),
        antimony: None,
    };
    let read = match sbml_io::read(xml) {
        Ok(r) => r,
        Err(diags) => {
            for d in diags {
                report.findings.push(Finding {
                    code: codes::SB_READ,
                    severity: Severity::Error,
                    element: Element::Model,
                    message: format!("Line {}, column {}: {}", d.line, d.column, d.message),
                    witness: None,
                    reaction: None,
                    title: None,
                    line: None,
                });
            }
            return report;
        }
    };
    let mut warnings = Vec::new();
    for d in read.warnings {
        warnings.push(Finding {
            code: codes::SB_READ_WARNING,
            severity: Severity::Info,
            element: Element::Model,
            message: format!("Line {}, column {}: {}", d.line, d.column, d.message),
            witness: None,
            reaction: None,
            title: None,
            line: None,
        });
    }
    let mut report = check_document(&read.document);
    report.findings.extend(warnings);
    sort(&mut report.findings);
    report
}

/// Runs every check on an SBML document already read (or built).
pub fn check_document(doc: &sbml_io::Document) -> Report {
    let mut report = Report {
        level: Some(doc.level),
        version: Some(doc.version),
        model_id: None,
        model_name: None,
        findings: Vec::new(),
        antimony: None,
    };
    for issue in sbml_io::validate(doc) {
        report.findings.push(Finding {
            code: codes::SB_STRUCTURE,
            severity: Severity::Error,
            element: Element::Model,
            message: issue.message,
            witness: None,
            reaction: None,
            title: None,
            line: None,
        });
    }
    let Some(model) = &doc.model else {
        report.findings.push(Finding {
            code: codes::SB_NO_MODEL,
            severity: Severity::Error,
            element: Element::Model,
            message: "The document contains no model.".into(),
            witness: None,
            reaction: None,
            title: None,
            line: None,
        });
        return report;
    };
    report.model_id = model.id.clone();
    report.model_name = model.name.clone();
    report.findings.extend(checks::run(model, doc.level, &doc.namespaces));
    report.antimony = Some(antimony_out::to_antimony(model, &doc.namespaces));
    // Show the reaction itself with every finding about one.
    for f in &mut report.findings {
        let id = match &f.element {
            Element::Reaction(id) => id,
            Element::LocalParameter { reaction, .. } => reaction,
            _ => continue,
        };
        if let Some(r) = model.reactions.iter().find(|r| &r.id == id) {
            f.reaction = Some(fmt::reaction(r, model));
        }
    }
    sort(&mut report.findings);
    report
}

/// Most severe first; within a severity, by code, then element, then line.
pub fn sort(findings: &mut [Finding]) {
    findings.sort_by(|a, b| {
        (a.severity, a.code, &a.element, a.line).cmp(&(b.severity, b.code, &b.element, b.line))
    });
}
