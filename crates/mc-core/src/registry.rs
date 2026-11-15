//! Database links (MIRIAM annotations) checked against the identifiers.org
//! registry ([`crate::registry_data`]): which database a link names, and
//! whether its id has that database's form.

use std::cell::RefCell;
use std::collections::HashMap;

use regex_lite::Regex;

use crate::registry_data::NAMESPACES;

/// One namespace of the registry.
#[derive(Debug, Clone, Copy)]
pub struct Namespace {
    pub prefix: &'static str,
    pub name: &'static str,
    pub pattern: &'static str,
    /// Its ids carry the prefix themselves: `CHEBI:15422`, `GO:0005829`.
    pub embedded: bool,
    pub deprecated: bool,
}

pub fn namespace(prefix: &str) -> Option<Namespace> {
    let p = prefix.to_ascii_lowercase();
    NAMESPACES
        .binary_search_by(|n| n.0.to_ascii_lowercase().cmp(&p))
        .ok()
        .map(|i| {
            let (prefix, name, pattern, embedded, deprecated) = NAMESPACES[i];
            Namespace {
                prefix,
                name,
                pattern,
                embedded,
                deprecated,
            }
        })
}

/// How a link is written.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Form {
    /// `https://identifiers.org/chebi/CHEBI:15422` or
    /// `https://identifiers.org/CHEBI:15422`.
    IdentifiersOrg,
    /// `urn:miriam:chebi:CHEBI%3A15422`: still understood, but deprecated.
    Urn,
    /// Any other address: it cannot be checked.
    Other,
}

/// A link taken apart.
#[derive(Debug, Clone, PartialEq)]
pub struct Link {
    pub form: Form,
    /// The database prefix as written, if the form has one.
    pub prefix: Option<String>,
    /// The id within the database, decoded.
    pub id: Option<String>,
}

/// `%3A` and the like decoded.
fn decode(s: &str) -> String {
    let bytes = s.as_bytes();
    let hex = |b: u8| (b as char).to_digit(16).map(|d| d as u8);
    let mut out = Vec::with_capacity(bytes.len());
    let mut i = 0;
    while i < bytes.len() {
        if bytes[i] == b'%'
            && i + 2 < bytes.len()
            && let (Some(h), Some(l)) = (hex(bytes[i + 1]), hex(bytes[i + 2]))
        {
            out.push(h * 16 + l);
            i += 3;
            continue;
        }
        out.push(bytes[i]);
        i += 1;
    }
    String::from_utf8_lossy(&out).into_owned()
}

pub fn parse(uri: &str) -> Link {
    let uri = uri.trim();
    if let Some(rest) = uri.strip_prefix("urn:miriam:") {
        let (prefix, id) = rest.split_once(':').unwrap_or((rest, ""));
        return Link {
            form: Form::Urn,
            prefix: Some(prefix.to_string()),
            id: (!id.is_empty()).then(|| decode(id)),
        };
    }
    let after = ["https://identifiers.org/", "http://identifiers.org/", "https://www.identifiers.org/", "http://www.identifiers.org/"]
        .iter()
        .find_map(|h| uri.strip_prefix(h));
    let Some(rest) = after else {
        return Link {
            form: Form::Other,
            prefix: None,
            id: None,
        };
    };
    let rest = decode(rest);
    // prefix/id: the older form, still everywhere.
    if let Some((prefix, id)) = rest.split_once('/') {
        return Link {
            form: Form::IdentifiersOrg,
            prefix: Some(prefix.to_string()),
            id: Some(id.to_string()),
        };
    }
    // Compact prefix:id. For a namespace whose ids carry the prefix
    // (CHEBI:15422), the whole is the id.
    let (prefix, id) = rest.split_once(':').unwrap_or((rest.as_str(), ""));
    let id = match namespace(prefix) {
        Some(ns) if ns.embedded => rest.clone(),
        _ => id.to_string(),
    };
    Link {
        form: Form::IdentifiersOrg,
        prefix: Some(prefix.to_string()),
        id: (!id.is_empty()).then_some(id),
    }
}

/// Compiled id patterns, kept for a run. A pattern the regex engine cannot
/// take (the registry's are written for Java) is not checked.
#[derive(Default)]
pub struct Patterns {
    cache: RefCell<HashMap<&'static str, Option<Regex>>>,
}

impl Patterns {
    /// Whether `id` has `ns`'s form; `None` if the pattern cannot be used.
    pub fn matches(&self, ns: &Namespace, id: &str) -> Option<bool> {
        let mut cache = self.cache.borrow_mut();
        let re = cache
            .entry(ns.prefix)
            .or_insert_with(|| Regex::new(ns.pattern).ok());
        re.as_ref().map(|r| r.is_match(id))
    }
}

/// What kind of thing a database's entries are, where that is clear.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    Molecule,
    Gene,
    Reaction,
    Literature,
    Taxon,
}

/// The kind of a well-known database's entries; `None` for databases whose
/// entries are of several kinds (GO, Reactome) or not listed here.
pub fn kind(prefix: &str) -> Option<Kind> {
    Some(match prefix.to_ascii_lowercase().as_str() {
        "chebi" | "kegg.compound" | "pubchem.compound" | "pubchem.substance" | "hmdb"
        | "kegg.drug" | "drugbank" | "bigg.metabolite" | "metanetx.chemical" | "lipidmaps"
        | "chembl.compound" | "inchikey" | "inchi" | "seed.compound" | "3dmet" => Kind::Molecule,
        "uniprot" | "ncbigene" | "ensembl" | "hgnc" | "kegg.genes" | "ncbiprotein" | "refseq"
        | "sgd" | "mgi" | "pr" | "interpro" | "pdb" | "wormbase" | "flybase" => Kind::Gene,
        "rhea" | "ec-code" | "kegg.reaction" | "bigg.reaction" | "metanetx.reaction"
        | "seed.reaction" | "sabiork.reaction" => Kind::Reaction,
        "pubmed" | "doi" | "arxiv" | "isbn" | "pmc" => Kind::Literature,
        "taxonomy" => Kind::Taxon,
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn links() {
        let l = parse("http://identifiers.org/chebi/CHEBI:15422");
        assert_eq!((l.form, l.prefix.as_deref(), l.id.as_deref()), (Form::IdentifiersOrg, Some("chebi"), Some("CHEBI:15422")));
        let l = parse("https://identifiers.org/CHEBI:15422");
        assert_eq!((l.prefix.as_deref(), l.id.as_deref()), (Some("CHEBI"), Some("CHEBI:15422")));
        let l = parse("https://identifiers.org/kegg.compound:C00002");
        assert_eq!((l.prefix.as_deref(), l.id.as_deref()), (Some("kegg.compound"), Some("C00002")));
        let l = parse("urn:miriam:obo.go:GO%3A0005829");
        assert_eq!((l.form, l.prefix.as_deref(), l.id.as_deref()), (Form::Urn, Some("obo.go"), Some("GO:0005829")));
        assert_eq!(parse("http://example.org/x").form, Form::Other);
    }

    #[test]
    fn patterns() {
        let p = Patterns::default();
        let chebi = namespace("chebi").unwrap();
        assert!(chebi.embedded);
        assert_eq!(p.matches(&chebi, "CHEBI:15422"), Some(true));
        assert_eq!(p.matches(&chebi, "15422"), Some(false));
        let uniprot = namespace("UniProt").unwrap();
        assert_eq!(p.matches(&uniprot, "P12345"), Some(true));
        assert!(namespace("refseq_synonym").is_none());
    }
}
