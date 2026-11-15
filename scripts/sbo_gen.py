# Generates crates/mc-core/src/sbo_data.rs from the SBO ontology:
# every term's id, name, parents and whether it is obsolete.
#
#   curl -LO https://raw.githubusercontent.com/EBI-BioModels/SBO/master/SBO_OBO.obo
#   python scripts/sbo_gen.py SBO_OBO.obo crates/mc-core/src/sbo_data.rs
import sys

path, out = sys.argv[1], sys.argv[2]
terms, cur, version = [], None, ""
for line in open(path, encoding="utf-8"):
    line = line.rstrip("\n")
    if line.startswith("property_value: owl:versionInfo"):
        version = line.split('"')[1]
    if line == "[Term]":
        cur = {"is_a": [], "obsolete": False}
        terms.append(cur)
        continue
    if line.startswith("["):
        cur = None
        continue
    if cur is None or ": " not in line:
        continue
    k, v = line.split(": ", 1)
    if k == "id":
        cur["id"] = int(v.split(":")[1])
    elif k == "name":
        cur["name"] = v
    elif k == "is_a":
        cur["is_a"].append(int(v.split(" ")[0].split(":")[1]))
    elif k == "is_obsolete":
        cur["obsolete"] = v.strip() == "true"
terms = [t for t in terms if "id" in t]
terms.sort(key=lambda t: t["id"])

BS = chr(92)
QUOTE = chr(34)


def rs(s):
    s = s.encode("ascii", "replace").decode("ascii")
    return QUOTE + s.replace(BS, BS + BS).replace(QUOTE, BS + QUOTE) + QUOTE


with open(out, "w", encoding="ascii", newline="\n") as f:
    f.write("//! The Systems Biology Ontology: every term's id, name, parents (is_a)\n")
    f.write("//! and whether it is obsolete. Generated from SBO_OBO.obo\n")
    f.write("//! (https://github.com/EBI-BioModels/SBO, master branch; the file's own\n")
    f.write("//! version string is " + version + ", but the file was last changed on\n")
    f.write("//! 2023-05-16, and EBI's Ontology Lookup Service lists the same terms)\n")
    f.write("//! by scripts/sbo_gen.py; do not edit by hand.\n\n")
    f.write("/// (id, name, parents, obsolete), sorted by id.\n")
    f.write("pub const TERMS: &[(u32, &str, &[u32], bool)] = &[\n")
    for t in terms:
        f.write("    (%d, %s, &[%s], %s),\n" % (
            t["id"], rs(t.get("name", "")),
            ", ".join(str(p) for p in t["is_a"]),
            "true" if t["obsolete"] else "false"))
    f.write("];\n")
print(len(terms), "terms, version", version)
