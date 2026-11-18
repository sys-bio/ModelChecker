//! NW008: a reaction that uses up one form of a common cofactor without
//! making another form of it, or makes one without using one up:
//! `Glucose + ATP -> G6P` has lost its ADP. Cofactors are recognised by
//! name (id or name, ignoring case, punctuation and a compartment suffix
//! such as `_c`). When the other side has a single placeholder species
//! (`X`), it is named as the likely missing form.
//!
//! Leaving a cofactor's partner out is sometimes done on purpose to keep a
//! model small, so this is information, or a warning when a placeholder
//! suggests the partner was meant to be there.

use std::collections::BTreeSet;

use super::Ctx;
use crate::codes::NW_COFACTOR;
use crate::finding::{Element, Severity};

/// One form of a cofactor: how it is shown, the names it goes by (as
/// [`normalise`] leaves them), and the forms it usually pairs with.
struct Form {
    shown: &'static str,
    names: &'static [&'static str],
    partners: &'static str,
}

const fn form(shown: &'static str, names: &'static [&'static str], partners: &'static str) -> Form {
    Form {
        shown,
        names,
        partners,
    }
}

/// The cofactors, each a family of forms that turn into one another.
const FAMILIES: &[&[Form]] = &[
    &[
        form("ATP", &["atp", "adenosinetriphosphate"], "ADP or AMP"),
        form("ADP", &["adp", "adenosinediphosphate"], "ATP"),
        form("AMP", &["amp", "adenosinemonophosphate"], "ATP or ADP"),
        form("cAMP", &["camp", "cyclicamp"], "ATP or AMP"),
    ],
    &[
        form("GTP", &["gtp", "guanosinetriphosphate"], "GDP or GMP"),
        form("GDP", &["gdp", "guanosinediphosphate"], "GTP"),
        form("GMP", &["gmp", "guanosinemonophosphate"], "GTP or GDP"),
        form("cGMP", &["cgmp", "cyclicgmp"], "GTP or GMP"),
    ],
    &[
        form("UTP", &["utp", "uridinetriphosphate"], "UDP or UMP"),
        form("UDP", &["udp", "uridinediphosphate"], "UTP"),
        form("UMP", &["ump", "uridinemonophosphate"], "UTP or UDP"),
    ],
    &[
        form("CTP", &["ctp", "cytidinetriphosphate"], "CDP or CMP"),
        form("CDP", &["cdp", "cytidinediphosphate"], "CTP"),
        form("CMP", &["cmp", "cytidinemonophosphate"], "CTP or CDP"),
    ],
    &[
        form(
            "NAD+",
            &["nad", "nadplus", "nadox", "nadoxidised", "nadoxidized"],
            "NADH",
        ),
        form("NADH", &["nadh", "nadred", "nadh2"], "NAD+"),
    ],
    &[
        form("NADP+", &["nadp", "nadpplus", "nadpox"], "NADPH"),
        form("NADPH", &["nadph", "nadpred", "nadph2"], "NADP+"),
    ],
    &[
        form("FAD", &["fad"], "FADH2"),
        form("FADH2", &["fadh2", "fadh"], "FAD"),
    ],
    &[
        form("FMN", &["fmn"], "FMNH2"),
        form("FMNH2", &["fmnh2", "fmnh"], "FMN"),
    ],
    &[
        form("GSH", &["gsh", "glutathione", "reducedglutathione"], "GSSG"),
        form(
            "GSSG",
            &[
                "gssg",
                "oxidisedglutathione",
                "oxidizedglutathione",
                "glutathionedisulfide",
                "glutathionedisulphide",
            ],
            "GSH",
        ),
    ],
];

/// Coenzyme A: free CoA pairs with any acyl-CoA (acetyl-CoA, succinyl-CoA,
/// ...), recognised by its name ending in "CoA".
const COA_FAMILY: usize = usize::MAX;
const COA_NAMES: &[&str] = &["coa", "coash", "coenzymea", "hscoa"];

/// Lower case, letters and digits only: "NAD+" is "nad", "acetyl-CoA" is
/// "acetylcoa".
fn normalise(s: &str) -> String {
    s.chars()
        .filter(char::is_ascii_alphanumeric)
        .map(|c| c.to_ascii_lowercase())
        .collect()
}

/// The ways a species might be named: its id and name, each also without a
/// short trailing segment that is probably a compartment (`ATP_c`,
/// `NADH_mito`).
fn candidates(id: &str, name: Option<&str>) -> Vec<String> {
    let mut out = Vec::new();
    for s in std::iter::once(id).chain(name) {
        out.push(normalise(s));
        if let Some((stem, last)) = s.rsplit_once('_')
            && !stem.is_empty()
            && (1..=4).contains(&last.len())
        {
            out.push(normalise(stem));
        }
    }
    out.retain(|c| !c.is_empty());
    out
}

/// A species recognised as a cofactor: its family, how to show it, and
/// what it usually pairs with.
struct Known {
    family: usize,
    shown: String,
    partners: &'static str,
}

fn identify(id: &str, name: Option<&str>) -> Option<Known> {
    let cands = candidates(id, name);
    for (i, family) in FAMILIES.iter().enumerate() {
        for f in *family {
            if cands.iter().any(|c| f.names.contains(&c.as_str())) {
                return Some(Known {
                    family: i,
                    shown: f.shown.to_string(),
                    partners: f.partners,
                });
            }
        }
    }
    if cands.iter().any(|c| COA_NAMES.contains(&c.as_str())) {
        return Some(Known {
            family: COA_FAMILY,
            shown: "CoA".to_string(),
            partners: "an acyl-CoA such as acetyl-CoA",
        });
    }
    if cands.iter().any(|c| c.len() > 3 && c.ends_with("coa")) {
        return Some(Known {
            family: COA_FAMILY,
            shown: name.unwrap_or(id).to_string(),
            partners: "CoA",
        });
    }
    None
}

/// A name that stands in for a species not yet decided: `X`, `Y2`, `unknown`.
/// Other single letters are not: `S + ATP -> P` names its product.
fn placeholder(id: &str, name: Option<&str>) -> bool {
    let n = normalise(name.filter(|n| !n.trim().is_empty()).unwrap_or(id));
    let stem = n.trim_end_matches(|c: char| c.is_ascii_digit());
    matches!(
        stem,
        "x" | "y" | "z" | "unknown" | "unk" | "missing" | "tbd"
    )
}

pub fn check(ctx: &mut Ctx) {
    let m = ctx.view.model;
    let name_of = |id: &str| {
        m.species
            .iter()
            .find(|s| s.id == id)
            .and_then(|s| s.name.as_deref())
    };
    for r in &m.reactions {
        let reactants: BTreeSet<&str> = r.reactants.iter().map(|x| x.species.as_str()).collect();
        let products: BTreeSet<&str> = r.products.iter().map(|x| x.species.as_str()).collect();
        // A cofactor alone on one side with nothing on the other is a
        // deliberate source or sink (`ATP ->`), not a forgotten partner.
        if reactants.is_empty() || products.is_empty() {
            continue;
        }
        let known = |side: &BTreeSet<&str>| -> Vec<Known> {
            side.iter()
                .filter_map(|id| identify(id, name_of(id)))
                .collect()
        };
        let used = known(&reactants);
        let made = known(&products);
        let families: BTreeSet<usize> = used.iter().chain(&made).map(|k| k.family).collect();
        for family in families {
            let in_used = used.iter().find(|k| k.family == family);
            let in_made = made.iter().find(|k| k.family == family);
            let (k, consumed, other_side) = match (in_used, in_made) {
                (Some(k), None) => (k, true, &products),
                (None, Some(k)) => (k, false, &reactants),
                _ => continue,
            };
            // A single placeholder on the other side is probably the
            // missing form.
            let stand_ins: Vec<&str> = other_side
                .iter()
                .copied()
                .filter(|id| identify(id, name_of(id)).is_none() && placeholder(id, name_of(id)))
                .collect();
            let guess = (stand_ins.len() == 1).then(|| stand_ins[0]);
            let (verb, lacks) = if consumed {
                ("uses up", "makes no")
            } else {
                ("makes", "uses up no")
            };
            let mut message = format!(
                "The reaction {verb} {} but {lacks} {}.",
                k.shown, k.partners
            );
            match guess {
                Some(x) => message.push_str(&format!(
                    " Is `{x}` meant to be {}? Naming it so would balance the cofactor.",
                    k.partners
                )),
                None => message.push_str(
                    " Often a forgotten species; fine if left out on purpose to keep the model \
                     small.",
                ),
            }
            let title = if consumed {
                format!("{} used up without {}", k.shown, k.partners)
            } else {
                format!("{} made without {}", k.shown, k.partners)
            };
            ctx.add_titled(
                NW_COFACTOR,
                if guess.is_some() {
                    Severity::Warning
                } else {
                    Severity::Info
                },
                Element::Reaction(r.id.clone()),
                title,
                message,
                None,
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn shown(id: &str) -> Option<String> {
        identify(id, None).map(|k| k.shown)
    }

    #[test]
    fn recognises_cofactors_by_name() {
        assert_eq!(shown("ATP").as_deref(), Some("ATP"));
        assert_eq!(shown("atp_c").as_deref(), Some("ATP"));
        assert_eq!(shown("NADH_mito").as_deref(), Some("NADH"));
        assert_eq!(shown("NADPH").as_deref(), Some("NADPH"));
        assert_eq!(shown("CoA").as_deref(), Some("CoA"));
        assert_eq!(shown("AcCoA").as_deref(), Some("AcCoA"));
        assert_eq!(
            identify("s1", Some("NAD+")).map(|k| k.shown).as_deref(),
            Some("NAD+")
        );
        assert_eq!(shown("G6P"), None);
        assert_eq!(shown("Glucose"), None);
    }

    #[test]
    fn placeholders() {
        assert!(placeholder("X", None));
        assert!(placeholder("X2", None));
        assert!(placeholder("unknown", None));
        assert!(!placeholder("G6P", None));
        assert!(!placeholder("P", None));
        assert!(!placeholder("Glucose", None));
    }
}
