//! Annotations: database links (MIRIAM) and model information.
//!
//! - MD001 species, reactions or compartments without database links;
//! - MD002 model information missing (creator, date, publication,
//!   organism, description);
//! - MD003 links that cannot work: an unknown database prefix, or an id
//!   not of the database's form (from the identifiers.org registry);
//! - MD004 old-style links (`urn:miriam:`) and deprecated databases;
//! - MD005 an `is` link to a database of the wrong kind for the element;
//! - MD006 elements with cryptic ids and no names;
//! - MD007 links that are not identifiers.org links, which cannot be
//!   checked.

use std::collections::BTreeMap;

use sbml_io::Meta;
use sbml_io::miriam::{Annotation, read_annotation};

use super::Ctx;
use crate::codes::*;
use crate::finding::{Element, Severity};
use crate::registry::{self, Form, Kind, Patterns};

/// When at least this many problems of one kind are found, they are listed
/// in one finding for the model.
const GROUP_AT: usize = 5;

struct Item<'m> {
    element: Element,
    /// "species", "reaction", ... for coverage.
    kind: &'static str,
    name: Option<&'m str>,
    id: &'m str,
    annotation: Annotation,
}

pub fn check(ctx: &mut Ctx, namespaces: &[(Option<String>, String)]) {
    let m = ctx.view.model;
    let read = |meta: &Meta| -> Annotation {
        meta.annotation
            .as_deref()
            .map(|a| read_annotation(a, namespaces))
            .unwrap_or_default()
    };
    let mut items: Vec<Item> = Vec::new();
    for c in &m.compartments {
        items.push(Item {
            element: Element::Compartment(c.id.clone()),
            kind: "compartment",
            name: c.name.as_deref(),
            id: &c.id,
            annotation: read(&c.meta),
        });
    }
    for s in &m.species {
        items.push(Item {
            element: Element::Species(s.id.clone()),
            kind: "species",
            name: s.name.as_deref(),
            id: &s.id,
            annotation: read(&s.meta),
        });
    }
    for r in &m.reactions {
        items.push(Item {
            element: Element::Reaction(r.id.clone()),
            kind: "reaction",
            name: r.name.as_deref(),
            id: &r.id,
            annotation: read(&r.meta),
        });
    }
    for p in &m.parameters {
        items.push(Item {
            element: Element::Parameter(p.id.clone()),
            kind: "parameter",
            name: p.name.as_deref(),
            id: &p.id,
            annotation: read(&p.meta),
        });
    }
    for f in &m.function_definitions {
        items.push(Item {
            element: Element::FunctionDefinition(f.id.clone()),
            kind: "function",
            name: f.name.as_deref(),
            id: &f.id,
            annotation: read(&f.meta),
        });
    }
    let model_annotation = read(&m.meta);

    coverage(ctx, &items);
    model_information(ctx, &model_annotation);
    links(ctx, &items, &model_annotation);
    names(ctx, &items);
}

fn plural(n: usize, one: &str, many: &str) -> String {
    format!("{n} {}", if n == 1 { one } else { many })
}

fn abbreviate(items: &[&str], n: usize) -> String {
    if items.len() <= n {
        items.join(", ")
    } else {
        format!("{}, ... ({} in all)", items[..n].join(", "), items.len())
    }
}

/// MD001.
fn coverage(ctx: &mut Ctx, items: &[Item]) {
    for (kind, one, many) in [
        ("species", "species", "species"),
        ("reaction", "reaction", "reactions"),
        ("compartment", "compartment", "compartments"),
    ] {
        let all: Vec<&Item> = items.iter().filter(|i| i.kind == kind).collect();
        let missing: Vec<&str> = all
            .iter()
            .filter(|i| {
                !i.annotation
                    .terms
                    .iter()
                    .any(|t| !t.model && !t.resources.is_empty())
            })
            .map(|i| i.id)
            .collect();
        if missing.is_empty() {
            continue;
        }
        let what = if all.len() == 1 {
            format!("The only {one} has no database links")
        } else if missing.len() == all.len() {
            format!(
                "None of the {} have database links",
                plural(all.len(), one, many)
            )
        } else {
            format!(
                "{} of {} have no database links",
                missing.len(),
                plural(all.len(), one, many)
            )
        };
        ctx.add(
            MD_NO_LINKS,
            Severity::Info,
            Element::Model,
            format!(
                "{what} (annotations saying what each one is, such as a ChEBI or UniProt entry): \
                 {}. Links make a model easier to understand, compare and reuse.",
                abbreviate(&missing, 15)
            ),
            None,
        );
    }
}

/// MD002.
fn model_information(ctx: &mut Ctx, a: &Annotation) {
    let m = ctx.view.model;
    let mut missing: Vec<&str> = Vec::new();
    if a.history.creators.is_empty() {
        missing.push("a creator");
    }
    if a.history.created.is_none() {
        missing.push("a creation date");
    }
    let publication = a.terms.iter().any(|t| {
        t.qualifier == "isDescribedBy"
            && t.resources.iter().any(|r| {
                let l = registry::parse(r);
                l.prefix.as_deref().and_then(registry::kind) == Some(Kind::Literature)
            })
    });
    if !publication {
        missing.push("a reference to its publication (PubMed or DOI)");
    }
    if !a.terms.iter().any(|t| t.qualifier == "hasTaxon") {
        missing.push("the organism it describes (taxon)");
    }
    if m.meta.notes.is_none() {
        missing.push("a description (notes)");
    }
    if missing.is_empty() {
        return;
    }
    let list = match missing.as_slice() {
        [a] => a.to_string(),
        _ => format!(
            "{} and {}",
            missing[..missing.len() - 1].join(", "),
            missing[missing.len() - 1]
        ),
    };
    ctx.add(
        MD_MODEL_INFO,
        Severity::Info,
        Element::Model,
        format!(
            "The model has no information on {list}. Repositories such as BioModels ask for these."
        ),
        None,
    );
}

/// MD003, MD004, MD005 and MD007.
fn links(ctx: &mut Ctx, items: &[Item], model: &Annotation) {
    let patterns = Patterns::default();
    let mut unknown: BTreeMap<String, usize> = BTreeMap::new();
    let mut malformed: Vec<(Element, String)> = Vec::new();
    let mut deprecated: BTreeMap<String, usize> = BTreeMap::new();
    let mut urns = 0usize;
    let mut other: Vec<String> = Vec::new();

    let all = items
        .iter()
        .map(|i| (&i.element, i.kind, &i.annotation))
        .chain(std::iter::once((&Element::Model, "model", model)));
    for (element, kind, annotation) in all {
        for term in &annotation.terms {
            for uri in &term.resources {
                let link = registry::parse(uri);
                match link.form {
                    Form::Other => {
                        other.push(uri.clone());
                        continue;
                    }
                    Form::Urn => urns += 1,
                    Form::IdentifiersOrg => {}
                }
                let Some(prefix) = link.prefix.as_deref() else {
                    continue;
                };
                // MIRIAM URNs named OBO ontologies `obo.go`, `obo.chebi`.
                let ns = registry::namespace(prefix)
                    .or_else(|| prefix.strip_prefix("obo.").and_then(registry::namespace));
                let Some(ns) = ns else {
                    *unknown.entry(prefix.to_string()).or_default() += 1;
                    continue;
                };
                if ns.deprecated {
                    *deprecated.entry(ns.prefix.to_string()).or_default() += 1;
                }
                // Common, resolvable variants are taken too: an ontology id
                // written with its own prefix (efo:EFO:0004541), and either
                // separator in a CURIE (the registry's OBI pattern wants
                // OBI_0000003; everyone writes OBI:0000003).
                let fits = |id: &str| {
                    let ok = |x: &str| patterns.matches(&ns, x) != Some(false);
                    ok(id)
                        || id
                            .split_once(':')
                            .filter(|(p, _)| p.eq_ignore_ascii_case(ns.prefix))
                            .is_some_and(|(_, rest)| ok(rest))
                        || ok(&id.replacen(':', "_", 1))
                        || ok(&id.replacen('_', ":", 1))
                };
                if let Some(id) = link.id.as_deref()
                    && !fits(id)
                {
                    malformed.push((
                        element.clone(),
                        format!(
                            "`{uri}`: `{id}` is not a {} id (the form is {})",
                            ns.name, ns.pattern
                        ),
                    ));
                }
                // MD005: `is` a thing of the wrong kind.
                if !term.model
                    && term.qualifier == "is"
                    && let Some(k) = registry::kind(prefix)
                {
                    let wrong = match kind {
                        "species" => matches!(k, Kind::Reaction | Kind::Literature | Kind::Taxon),
                        "reaction" => matches!(
                            k,
                            Kind::Molecule | Kind::Gene | Kind::Literature | Kind::Taxon
                        ),
                        "compartment" => {
                            matches!(k, Kind::Molecule | Kind::Reaction | Kind::Literature)
                        }
                        _ => false,
                    };
                    if wrong {
                        let what = match k {
                            Kind::Molecule => "a molecule",
                            Kind::Gene => "a gene or protein",
                            Kind::Reaction => "a reaction",
                            Kind::Literature => "a publication",
                            Kind::Taxon => "an organism",
                        };
                        ctx.add(
                            MD_WRONG_KIND,
                            Severity::Warning,
                            element.clone(),
                            format!(
                                "The {kind} is annotated as being (`is`) {what}: `{uri}` ({}). \
                                 A {kind} cannot be {what}; perhaps another qualifier was meant \
                                 (such as `hasPart`, `isVersionOf` or `isDescribedBy`), or the \
                                 link belongs on another element.",
                                ns.name
                            ),
                            None,
                        );
                    }
                }
            }
        }
    }

    if !unknown.is_empty() {
        let total: usize = unknown.values().sum();
        let list = unknown
            .iter()
            .map(|(p, n)| format!("{p} ({n})"))
            .collect::<Vec<_>>()
            .join(", ");
        ctx.add(
            MD_BROKEN_LINK,
            Severity::Warning,
            Element::Model,
            format!(
                "{} use{} a database prefix that identifiers.org does not know, so {} cannot be \
                 followed: {list}.",
                plural(total, "link", "links"),
                if total == 1 { "s" } else { "" },
                if total == 1 { "it" } else { "they" }
            ),
            None,
        );
    }
    if malformed.len() >= GROUP_AT {
        let examples: Vec<String> = malformed
            .iter()
            .take(8)
            .map(|(e, s)| format!("{e}: {s}"))
            .collect();
        ctx.add(
            MD_BROKEN_LINK,
            Severity::Warning,
            Element::Model,
            format!(
                "{} have an id that is not in their database's form, so they cannot be followed. \
                 For example:\n{}",
                plural(malformed.len(), "link", "links"),
                examples.join("\n")
            ),
            None,
        );
    } else {
        for (e, s) in malformed {
            ctx.add(
                MD_BROKEN_LINK,
                Severity::Warning,
                e,
                format!("A database link has an id that is not in its database's form: {s}."),
                None,
            );
        }
    }
    if urns > 0 || !deprecated.is_empty() {
        let mut parts = Vec::new();
        if urns > 0 {
            parts.push(format!(
                "{} use{} the old `urn:miriam:` form, which identifiers.org has replaced with \
                 `https://identifiers.org/prefix:id`",
                plural(urns, "link", "links"),
                if urns == 1 { "s" } else { "" }
            ));
        }
        if !deprecated.is_empty() {
            let list = deprecated
                .iter()
                .map(|(p, n)| format!("{p} ({n})"))
                .collect::<Vec<_>>()
                .join(", ");
            parts.push(format!(
                "links point to databases identifiers.org lists as deprecated: {list}"
            ));
        }
        ctx.add(
            MD_OLD_LINK,
            Severity::Info,
            Element::Model,
            format!("{}. They still work for now.", parts.join("; ")),
            None,
        );
    }
    if !other.is_empty() {
        let shown: Vec<&str> = other.iter().map(String::as_str).collect();
        ctx.add(
            MD_OTHER_LINK,
            Severity::Info,
            Element::Model,
            format!(
                "{} {} not identifiers.org links, so {} cannot be checked: {}.",
                plural(other.len(), "link", "links"),
                if other.len() == 1 { "is" } else { "are" },
                if other.len() == 1 { "it" } else { "they" },
                abbreviate(&shown, 5)
            ),
            None,
        );
    }
}

/// An id that tells a reader nothing: `species_0000012`, `s12`, `re5`, or
/// a long machine-made one.
fn cryptic(id: &str) -> bool {
    let digits_tail = id.trim_end_matches(|c: char| c.is_ascii_digit());
    let stem = digits_tail.trim_end_matches('_').to_ascii_lowercase();
    // Generic words, and CellDesigner's s12, sa12, re12. Short hand-made
    // ids such as J1, R1 or v1 read fine in context and are not flagged.
    let generic = matches!(
        stem.as_str(),
        "s" | "sa" | "re" | "species" | "reaction" | "compartment" | "metabolite" | "parameter"
    );
    let numbered = digits_tail.len() < id.len();
    let machine =
        id.len() > 24 && id.chars().filter(|c| c.is_ascii_hexdigit()).count() * 3 > id.len() * 2;
    (generic && numbered) || machine
}

/// MD006.
fn names(ctx: &mut Ctx, items: &[Item]) {
    let unnamed: Vec<&str> = items
        .iter()
        .filter(|i| matches!(i.kind, "species" | "reaction" | "compartment"))
        .filter(|i| i.name.is_none_or(|n| n.trim().is_empty()) && cryptic(i.id))
        .map(|i| i.id)
        .collect();
    if unnamed.is_empty() {
        return;
    }
    ctx.add(
        MD_NO_NAME,
        Severity::Info,
        Element::Model,
        format!(
            "{} {} only a machine-style id and no name, which makes the model hard to read: {}.",
            plural(unnamed.len(), "element", "elements"),
            if unnamed.len() == 1 { "has" } else { "have" },
            abbreviate(&unnamed, 12)
        ),
        None,
    );
}

#[cfg(test)]
mod tests {
    use super::cryptic;

    #[test]
    fn cryptic_ids() {
        assert!(cryptic("species_0000012"));
        assert!(cryptic("s12"));
        assert!(cryptic("re5"));
        assert!(cryptic("mwf626e95e_543f_41e4_aad4_c6bf60ab345b"));
        assert!(!cryptic("ATP"));
        assert!(!cryptic("glucose"));
        assert!(!cryptic("J1_hexokinase"));
    }
}
