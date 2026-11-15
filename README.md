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
- **Context.** Each finding about a reaction shows the reaction itself,
  e.g. `J4: S4 -> ; k4*(S4 + E)`, with `$` marking boundary species.
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
