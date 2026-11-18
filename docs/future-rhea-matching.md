# Possible future feature: matching reactions to Rhea

Status: idea, not started. Written November 2026, after version 1.0.2 added
NW008 (a cofactor without its partner).

## The idea

Rhea (https://www.rhea-db.org) is an expert-curated database of biochemical
reactions. Every participant is a ChEBI compound, and every reaction is
balanced for mass and charge. If the checker could match a model's
reactions to Rhea, it could:

- **Fill in missing species.** `Glucose + ATP -> G6P + X` matches the
  hexokinase reaction (D-glucose + ATP = D-glucose 6-phosphate + ADP + H+),
  so `X` is ADP. NW008 already catches this for common cofactors from a
  built-in list. Rhea would cover any participant, such as a forgotten
  phosphate, CO2 or a second product.
- **Spot wrong or extra species.** A reaction close to a Rhea reaction but
  differing by one species is probably a typo or a lumping choice worth a
  note.
- **Check stoichiometry and direction** against the curated reaction.
- **Name the reaction.** It could suggest the enzyme (EC number) and add a
  Rhea link (`https://identifiers.org/rhea:NNNNN`) as an annotation, which
  ties in with the MD (annotation) checks.

## How it could work

1. **Identify each species as a ChEBI compound.**
   - Exact: from the species' own ChEBI annotation, when it has one.
   - Approximate: from its id or name against ChEBI names and synonyms
     ("Glucose", "glc", "D-glucose"), after the same normalising NW008 does
     (case, punctuation, compartment suffixes such as `_c`).
   - Different charge states of one compound (ATP and ATP(4-)) count as the
     same. ChEBI records these as conjugate acids and bases.
2. **Match each reaction** as a set of ChEBI compounds against Rhea,
   ignoring H+ and H2O, which kinetic models usually leave out.
   - Full match: the reaction is a known one (information, with the link).
   - Match except for one placeholder species: that species is probably the
     remaining Rhea participant (warning, with the suggestion).
   - Near match (one species extra or missing): a warning naming the
     difference.
3. **Report** through a new finding group, for example:
   - RH001 matches Rhea reaction N (information, with link and EC number);
   - RH002 close to Rhea reaction N but differs by these species (warning);
   - RH003 placeholder `X` is probably compound Y (warning).

## Keeping the model on the user's computer

The page promises that the model never leaves the user's computer. Three
ways to get the data, and how each keeps that promise:

| Option | How | Privacy | Cost |
|---|---|---|---|
| A. Built into the page | A compact index of Rhea and the ChEBI names it uses, compiled into the WebAssembly | Kept | Larger download for every user, even those who never use it |
| B. Separate data file (recommended) | The same index as a file published with the site, fetched only when the user asks for Rhea matching, then cached by the browser | Kept: only the data file is downloaded, nothing is sent | One extra download, only for users who want it |
| C. Online queries | Ask Rhea's web service about each reaction | Broken: species names go to EBI | No data to ship, but needs the internet and an opt-in |

Option B fits best. The index would be made by a script in `scripts/`, like
`registry_gen.py`, from Rhea's and ChEBI's bulk downloads, and rebuilt when
Rhea is updated.

## Difficulties to expect

- **Names are ambiguous.** "Glucose" could be D-glucose, alpha-D-glucose or
  a generic hexose, and "G6P" is an abbreviation ChEBI may not list. Name
  matching will need a small hand-made table of common abbreviations (G6P,
  F6P, PEP, DHAP, GAP, ...) on top of the ChEBI synonyms.
- **Kinetic models lump reactions.** One model reaction often stands for
  several steps (glycolysis in three reactions), which no single Rhea
  reaction matches. These should get no finding, not a false warning.
- **Abstract models** (`S1 -> S2`) have nothing to match and must be left
  alone.
- **Generic compounds.** Rhea has reactions on classes ("an aldose"), which
  can match too broadly.
- **Transport reactions** (the same compound in two compartments) need
  their own handling.
- **Size.** The index size has to be measured. Restricting ChEBI to
  compounds that appear in Rhea, and storing ids compactly, should keep it
  reasonable, but this is the first thing to test.

## Licensing

Rhea and ChEBI data are distributed under the Creative Commons Attribution
licence (CC BY 4.0). This should be confirmed when the work starts. The page
and the README would need to credit them.

## First steps, when this is picked up

1. Build the index from the Rhea and ChEBI downloads and measure its size.
2. Measure on a collection of curated models (for example BioModels): how
   many species can be identified, how many reactions match, and how often
   the suggestions are right.
3. Only then decide on the finding codes and the page's controls.
