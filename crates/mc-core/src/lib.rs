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

pub mod checks;
pub mod codes;
pub mod eval;
pub mod finding;
pub mod fmt;
pub mod json;
pub mod moiety;
pub mod sample;
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
                });
            }
            return report;
        }
    };
    let doc = read.document;
    report.level = Some(doc.level);
    report.version = Some(doc.version);
    for d in read.warnings {
        report.findings.push(Finding {
            code: codes::SB_READ_WARNING,
            severity: Severity::Info,
            element: Element::Model,
            message: format!("Line {}, column {}: {}", d.line, d.column, d.message),
            witness: None,
        });
    }
    for issue in sbml_io::validate(&doc) {
        report.findings.push(Finding {
            code: codes::SB_STRUCTURE,
            severity: Severity::Error,
            element: Element::Model,
            message: issue.message,
            witness: None,
        });
    }
    let Some(model) = &doc.model else {
        report.findings.push(Finding {
            code: codes::SB_NO_MODEL,
            severity: Severity::Error,
            element: Element::Model,
            message: "The document contains no model.".into(),
            witness: None,
        });
        return report;
    };
    report.model_id = model.id.clone();
    report.model_name = model.name.clone();
    report.findings.extend(checks::run(model, doc.level));
    // Most severe first; within a severity, by code, then element.
    report
        .findings
        .sort_by(|a, b| (a.severity, a.code, &a.element).cmp(&(b.severity, b.code, &b.element)));
    report
}
