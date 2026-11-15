# Generates crates/mc-core/src/registry_data.rs from the identifiers.org
# registry: every namespace's prefix, name, id pattern, whether its ids
# carry their own prefix (CHEBI:15422), and whether it is deprecated.
#
#   curl -L -o resolver.json https://registry.api.identifiers.org/resolutionApi/getResolverDataset
#   python scripts/registry_gen.py resolver.json crates/mc-core/src/registry_data.rs
import json
import sys
from datetime import date

src, out = sys.argv[1], sys.argv[2]
data = json.load(open(src, encoding="utf-8"))
namespaces = sorted(data["payload"]["namespaces"], key=lambda n: n["prefix"].lower())

QUOTE = chr(34)
BS = chr(92)


def rs(s):
    s = (s or "").encode("ascii", "replace").decode("ascii")
    return QUOTE + s.replace(BS, BS + BS).replace(QUOTE, BS + QUOTE) + QUOTE


with open(out, "w", encoding="ascii", newline="\n") as f:
    f.write("//! The identifiers.org registry: every namespace's prefix, name, id\n")
    f.write("//! pattern, whether its ids carry their own prefix (CHEBI:15422), and\n")
    f.write("//! whether it is deprecated. Generated from the registry's resolver\n")
    f.write("//! dataset (https://registry.api.identifiers.org/resolutionApi/getResolverDataset,\n")
    f.write("//! downloaded " + date.today().isoformat() + ") by scripts/registry_gen.py; do not edit by hand.\n\n")
    f.write("/// (prefix, name, pattern, prefix embedded in ids, deprecated), sorted\n")
    f.write("/// by lower-case prefix.\n")
    f.write("pub const NAMESPACES: &[(&str, &str, &str, bool, bool)] = &[\n")
    for n in namespaces:
        f.write("    (%s, %s, %s, %s, %s),\n" % (
            rs(n["prefix"]), rs(n["name"]), rs(n["pattern"]),
            "true" if n.get("namespaceEmbeddedInLui") else "false",
            "true" if n.get("deprecated") else "false"))
    f.write("];\n")
print(len(namespaces), "namespaces")
