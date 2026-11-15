# Model Checker

## AI Assistance

This software was developed with the assistance of Anthropic's Claude, an AI model, through Claude Code. Stop here if you aren't comfortable with AI assisted code generation.

## About

Static checks for kinetic SBML models. It finds rate laws that divide by zero,
reactions that can drive species negative, zero-order reactions that are
probably typos, and undefined or unused symbols. It reports each problem as a
finding, with a concrete witness where one exists.

**Try it: https://sys-bio.github.io/ModelChecker/**

It runs in the browser (Rust compiled to WebAssembly): the model never leaves
your computer. There is also a command-line tool for checking whole folders.

## What is checked

| Code  | Severity | Check |
|-------|----------|-------|
| RL012 | error    | Rate law or rule singular at the initial values (simulation fails at once) |
| RL001 | warning  | Rate law singular at some state with all species >= 0 |
| RL002 | info     | Possibly singular: interval analysis could not rule it out |
| RL003 | warning  | Species can go negative: with the species at zero, the reaction still consumes it |
| RL004 | info     | Species may go negative: not proved, or fails only outside initial values |
| RL005 | warning  | Zero order: the reaction consumes species its rate does not depend on |
| RL006 | warning  | A reactant is missing from the rate law |
| RL007 | info     | Constant-rate source |
| RL008 | warning  | No kinetic law |
| RL009 | warning  | Species used in a rate law but not declared as reactant, product or modifier |
| RL010 | info     | Modifier not used by the rate law |
| RL011 | warning  | Reaction with no reactants and no products |
| RL013 | warning  | Mass-action rate law whose powers differ from the stoichiometry (e.g. `k*A` for `2 A -> B`) |
| RU001 | error    | Assignment rules depending on each other in a loop (in Level 3, possibly through a reaction's rate) |
| RU002 | error    | Initial assignments (with rules) depending on each other in a loop |
| RU003 | error    | Variable set more than once: two rules, a rule and an initial assignment, a rule and an event, or reactions and a rule on a non-boundary species |
| RU004 | error    | Something marked constant that a rule or event sets |
| RU005 | info     | Algebraic rules (many simulators do not support them) |
| RU006 | warning  | Event trigger already true at the initial values (says whether the event then fires at time 0 or waits) |
| RU007 | warning  | Event that can never fire: its trigger depends on nothing that changes |
| RU008 | info     | Event with no assignments |
| RU009 | warning  | Negative event delay |
| NW001 | info     | Species only produced, never consumed (accumulates) |
| NW002 | info     | Species only consumed, never produced (runs down) |
| NW003 | info     | Species in no reaction, rule or event, and used nowhere |
| NW004 | info     | Network in separate parts that never affect each other |
| NW005 | info     | Conserved totals (conservation laws), e.g. `MKKK + MKKK_P = 100` |
| NW006 | warning  | Species that start at zero and can never be made (an empty siphon), and the reactions that therefore never run |
| NW007 | info     | Reaction that changes no variable species (only carries a flux) |
| TD001 | warning  | Reversible rate law whose equilibrium is not a constant mass-action ratio (Haldane), e.g. shifted by a modifier |
| TD002 | warning  | Cycle of reversible reactions whose equilibrium constants do not multiply to 1 (Wegscheider) |
| UN001 | warning  | Rate law units differ from extent per time |
| UN002 | warning  | Adds, subtracts or compares quantities of different dimensions |
| UN003 | warning  | Same, with the same dimensions but a different scale (mM vs M) |
| UN004 | warning  | exp, ln, trig or an exponent applied to a quantity with units |
| UN005 | warning  | Rate law in concentration per time: missing a compartment volume factor |
| UN006 | warning  | Same suspicion without units: concentrations in a compartment of size other than 1, and no compartment in the rate law |
| UN007 | warning  | Rule or assignment units differ from its variable's |
| UN008 | info     | What could not be checked (usually parameters without units) |
| SO001 | warning  | SBO term that is not valid or not in the ontology |
| SO002 | info     | Obsolete SBO term |
| SO003 | warning  | SBO term from the wrong branch for the element (e.g. a reaction term on a kinetic law) |
| SO004 | warning  | Kinetic law marked irreversible, but the rate can be negative |
| SO005 | warning  | Kinetic law marked reversible, but the rate can never be negative |
| SO006 | warning  | Kinetic law marked mass action, but not of mass-action form |
| SO007 | warning  | Mass-action term whose orders or reactant/product counts differ from the rate law or the reaction |
| SO008 | warning  | Enzymatic or Hill term that does not fit: no saturation, no species power, no modifier, wrong reactant count |
| AN001 | error    | The edited Antimony text cannot be read (web page only) |
| AN002 | warning  | Antimony warning on the edited text, e.g. a name never declared (web page only) |
| SY001 | error    | Undefined symbol |
| SY002-SY004 | info | Unused parameter, local parameter, function |
| SY005 | info     | Local parameter shadows a global |
| SY006 | warning  | Symbol with no value from anywhere |
| SY007 | info     | MathML that could not be analysed |
| SB001-SB004 | | Reading problems, structural problems (duplicate ids, dangling references) |

### How

- Function definitions and assignment rules are expanded into the rate laws
  that use them, so dependencies through them count.
- Species are non-negative (unless they start negative, like a voltage), and
  parameters are fixed at their values unless a rule or event changes them.
- **Singularities.** Each expression is evaluated over intervals, which
  covers all such states at once. A concrete point is then searched for. A
  hazard inside a `piecewise` branch whose condition may guard it is not
  reported.
- **Positivity.** For each species S with net stoichiometry n, the rate v
  must satisfy n * v >= 0 when S = 0. A reversible rate law that runs
  backwards at S = 0 (and so produces S) passes. The `reversible` attribute
  is ignored: only the rate law decides which way a reaction can run. A
  counterexample is first sought with the other variables between zero and
  their initial values. A failure found only outside that range is reported
  as info, since it is often a part exceeding its total (`ET - ES`).
- **Network.** Conservation laws come from the left null space of the
  stoichiometry matrix, in amounts, over the species only reactions change.
  The basis is then recombined so the laws read naturally: as few negative
  coefficients and species as possible, so `S + ES + P` rather than
  `-E + S + P`. Species that can never appear are found by a fixed point.
  Start with every species at zero, then free any species that a reaction
  can make while the remaining ones are still zero. Whether a reaction can
  run, and which way, is asked of its rate law (by interval arithmetic), so a
  rate law that ignores its reactant still runs. Whatever is left can never
  become non-zero. When many species share a network finding, they are
  listed in one finding for the model.
- **Thermodynamics.** Many published models are not thermodynamically
  consistent, often harmlessly, so these are warnings. For each reversible
  rate law that depends on both its reactants and its products, the
  mass-action ratio at which it is zero is found numerically at several
  random states, by scaling the products until the rate vanishes. A
  consistent rate law gives one value, its equilibrium constant. If the value
  moves, each variable is perturbed in turn to find what moves it (TD001).
  Cycles of reversible reactions come from the null space of the
  stoichiometry, counting boundary species, so a path between two boundary
  species is not a cycle. Their equilibrium constants must multiply to 1
  around each cycle (TD002). Products within 1% of 1 are taken as rounding.
  Only reactions with whole-number stoichiometries up to 4 are considered.
- **Units.** Each expression's units are inferred from the declared units
  of its symbols (SI exponents plus a scale). A species is substance over
  its compartment's size unless `hasOnlySubstanceUnits`. A parameter without
  units is unknown, unless an assignment rule sets it, in which case it takes
  the rule's units. Unknowns make a product unknown, so unit findings arise
  only where the model declares enough. A bare number takes the units of
  whatever it is added to. In a product, only small whole numbers and 1/2
  count as dimensionless; other numbers may hide units (conversion factors,
  rate constants), as in libSBML.
- **SBO terms.** The ontology is built in (`crates/mc-core/src/sbo_data.rs`,
  generated from the EBI-BioModels SBO repository by
  `scripts/sbo_gen.py`). Each element's term must come from the branch SBML
  allows for it. A kinetic law's term is checked against its formula, as
  written with function definitions inlined. Mass-action terms are checked by
  writing the rate law as signed products of species powers and comparing the
  forward and reverse parts with the orders and counts in the term's name and
  with the reaction's stoichiometry. Declared modifiers, such as a template or
  catalyst, may appear without counting as a stoichiometry mismatch.
- **Repeated findings.** When three or more rate laws fail a unit check for
  the same reason (off by the same ratio, or missing the same compartment's
  volume), they are reported as one model-level finding that lists the
  reactions and suggests the shared cause.
- **Mass-action order.** A rate law of mass-action form (one forward term,
  at most one reverse term) should raise each species to its stoichiometry.
  Only plausibly elementary reactions are checked: whole-number
  stoichiometries, at most three molecules a side, and no power-law
  (GMA/S-system) exponents.
- **Context.** Each finding about a reaction shows the reaction itself,
  e.g. `J4: S4 -> ; k4*(S4 + E)`, with `$` marking boundary species. The
  page also shows the whole model as Antimony beside the findings. Click a
  finding to see it in the model, or click a marked line to see its
  findings. The panel is editable. Change the model and press **Re-check**
  (or Ctrl+Enter) to check the edited version, with findings on your lines.
  Antimony errors and warnings, such as a name used but never declared, are
  reported too (AN001, AN002). **Revert to original** returns to the loaded
  file, and **Download SBML** saves the edited model.
- **Editing never loses anything.** The loaded SBML stays the master copy.
  An edit of the text is merged into it by comparing the edited text with
  the text as first written: only what you changed is changed. The text
  shows unit definitions (`unit mM = 0.001 mole litre^-1`), units of
  variables (`k1 has per_s`) and SBO terms (`J1.sboTerm = 176`), and these
  can be edited. Model-wide units, rate-law SBO terms, and units whose id is
  an Antimony keyword have no Antimony syntax, so they are shown as comments.
  Annotations, notes and names are not shown. Everything the text does not
  carry is kept, as are untouched elements, byte for byte. Local parameters
  are written as `J1_k1` and turned back into local parameters, and numbers
  are written exactly. Tested on all BioModels
  (`cargo run --release -p mc-core --features antimony --example edit_roundtrip -- FOLDER`):
  an unedited text gives back the original exactly (1075 of 1075). Changing
  one value, one rate law, one SBO term or one unit changes only that (928 of
  928, 953 of 953, 335 of 335, 328 of 328). Ids that are Antimony keywords, such as a compartment called
  `compartment`, are written with a trailing `_`. Every BioModels model
  converts to Antimony that the antimony crate parses
  (`cargo run --release -p mc-core --example ant_roundtrip -- FOLDER`).
- **Conservation.** A witness that empties a whole conserved moiety (found
  from the stoichiometry matrix) is not a reachable state, so it is not used.

## Layout

```
crates/
  sbml-io/    SBML reader (shared with Annotator and PathwayDesigner)
  antimony/   Antimony to SBML, used by tests and the ant2sbml example
  mc-core/    the checks: no I/O, builds for wasm32
  mc-cli/     `modelcheck`, the command-line tool
  mc-web/     the web page (Trunk + wasm-bindgen)
```

## Building

```
cargo test -p mc-core                      # unit and model tests
cargo run --release -p mc-cli -- model.xml # check one model
cargo run --release -p mc-cli -- --summary folder/   # counts over a corpus
cargo run -p mc-core --example ant2sbml -- model.ant > model.xml
```

The web page:

```
cd crates/mc-web
trunk build --release --public-url ./      # output in crates/mc-web/dist
trunk serve                                # local preview
```

If `trunk` fails with `invalid value '1' for '--no-color'`, set
`NO_COLOR=true` first.

Pushing to `main` builds and deploys the page to GitHub Pages
(`.github/workflows/ci.yml`). This needs a one-time setting on GitHub:
Settings -> Pages -> Source: GitHub Actions.

## Funding

This work was supported by NIH Biomedical Imaging and Bioengineering award
P41 EB023912 through HMS at the Center for Reproducible Biomedical Modeling,
[reproduciblebiomodels.org](https://reproduciblebiomodels.org).
