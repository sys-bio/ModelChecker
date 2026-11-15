//! What a check reports.

use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Severity {
    /// The model is wrong: a concrete point was found where it fails.
    Error,
    /// Probably wrong, or could not be shown to be right.
    Warning,
    /// Worth knowing; often intended.
    Info,
}

impl Severity {
    pub fn as_str(self) -> &'static str {
        match self {
            Severity::Error => "error",
            Severity::Warning => "warning",
            Severity::Info => "info",
        }
    }
}

impl fmt::Display for Severity {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

/// The model element a finding is about.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Element {
    Model,
    Reaction(String),
    Species(String),
    Parameter(String),
    LocalParameter { reaction: String, id: String },
    Compartment(String),
    FunctionDefinition(String),
    /// An assignment or rate rule, named by the variable it sets.
    Rule(String),
    /// An algebraic rule, by position.
    AlgebraicRule(usize),
    InitialAssignment(String),
    /// An event, by id if it has one, else by position.
    Event(String),
}

impl Element {
    pub fn kind(&self) -> &'static str {
        match self {
            Element::Model => "model",
            Element::Reaction(_) => "reaction",
            Element::Species(_) => "species",
            Element::Parameter(_) => "parameter",
            Element::LocalParameter { .. } => "local parameter",
            Element::Compartment(_) => "compartment",
            Element::FunctionDefinition(_) => "function",
            Element::Rule(_) => "rule",
            Element::AlgebraicRule(_) => "algebraic rule",
            Element::InitialAssignment(_) => "initial assignment",
            Element::Event(_) => "event",
        }
    }

    pub fn id(&self) -> String {
        match self {
            Element::Model => String::new(),
            Element::Reaction(s)
            | Element::Species(s)
            | Element::Parameter(s)
            | Element::Compartment(s)
            | Element::FunctionDefinition(s)
            | Element::Rule(s)
            | Element::InitialAssignment(s)
            | Element::Event(s) => s.clone(),
            Element::LocalParameter { reaction, id } => format!("{reaction}.{id}"),
            Element::AlgebraicRule(i) => format!("#{}", i + 1),
        }
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Element::Model => f.write_str("model"),
            _ => write!(f, "{} {}", self.kind(), self.id()),
        }
    }
}

/// One problem found in a model.
#[derive(Debug, Clone, PartialEq)]
pub struct Finding {
    /// Stable code, for filtering and suppression. See [`crate::codes`].
    pub code: &'static str,
    pub severity: Severity,
    pub element: Element,
    pub message: String,
    /// A concrete point that shows the problem, for example
    /// `S1 = 0, P = 0, E = 1 gives v = 2.3`.
    pub witness: Option<String>,
    /// The reaction the finding is about, written out:
    /// `J1: S1 + 2 S2 -> P; k1*S1*S2`.
    pub reaction: Option<String>,
}

impl fmt::Display for Finding {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} {} [{}]: {}",
            self.severity, self.code, self.element, self.message
        )?;
        if let Some(w) = &self.witness {
            write!(f, " (witness: {w})")?;
        }
        Ok(())
    }
}
