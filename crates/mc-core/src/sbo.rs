//! Looking up Systems Biology Ontology terms: names, obsolescence, and
//! the is_a hierarchy. The data is in [`crate::sbo_data`].

use crate::sbo_data::TERMS;

/// Well-known terms.
pub const RATE_LAW: u32 = 1;
pub const PARTICIPANT_ROLE: u32 = 3;
pub const MODELLING_FRAMEWORK: u32 = 4;
pub const MASS_ACTION: u32 = 12;
pub const MODIFIER: u32 = 19;
pub const MATHEMATICAL_EXPRESSION: u32 = 64;
pub const HILL: u32 = 192;
pub const OCCURRING_ENTITY: u32 = 231;
pub const PHYSICAL_ENTITY: u32 = 236;
pub const SYSTEMS_DESCRIPTION_PARAMETER: u32 = 545;
/// Irreversible non-modulated non-interacting unireactant enzymes: the
/// Michaelis-Menten family.
pub const MICHAELIS_MENTEN_FAMILY: u32 = 28;
pub const ENZYMATIC_MODULATED: u32 = 430;
pub const REVERSIBLE_HILL: u32 = 425;
pub const MODULAR: u32 = 527;

#[derive(Debug, Clone, Copy)]
pub struct Term {
    pub id: u32,
    pub name: &'static str,
    pub parents: &'static [u32],
    pub obsolete: bool,
}

/// `SBO:0000049` (or `0000049`, or `49`) to 49.
pub fn parse(s: &str) -> Option<u32> {
    let s = s.trim();
    let digits = s.strip_prefix("SBO:").unwrap_or(s);
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    digits.parse().ok()
}

pub fn format(id: u32) -> String {
    format!("SBO:{id:07}")
}

pub fn term(id: u32) -> Option<Term> {
    TERMS
        .binary_search_by_key(&id, |t| t.0)
        .ok()
        .map(|i| {
            let (id, name, parents, obsolete) = TERMS[i];
            Term {
                id,
                name,
                parents,
                obsolete,
            }
        })
}

/// Whether `id` is `ancestor` or below it.
pub fn is_a(id: u32, ancestor: u32) -> bool {
    let mut stack = vec![id];
    let mut seen = Vec::new();
    while let Some(t) = stack.pop() {
        if t == ancestor {
            return true;
        }
        if seen.contains(&t) {
            continue;
        }
        seen.push(t);
        if let Some(term) = term(t) {
            stack.extend_from_slice(term.parents);
        }
    }
    false
}

/// `SBO:0000049 (mass action rate law for ...)`.
pub fn describe(id: u32) -> String {
    match term(id) {
        Some(t) => format!("{} ({})", format(id), t.name),
        None => format(id),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hierarchy() {
        assert_eq!(parse("SBO:0000049"), Some(49));
        assert_eq!(parse("SBO:x"), None);
        assert!(is_a(49, RATE_LAW));
        assert!(is_a(49, MASS_ACTION));
        assert!(is_a(29, MICHAELIS_MENTEN_FAMILY));
        assert!(!is_a(176, RATE_LAW));
        assert_eq!(term(1).unwrap().name, "rate law");
        assert!(is_a(RATE_LAW, MATHEMATICAL_EXPRESSION));
    }
}
