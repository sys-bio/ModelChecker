# Model Checker

Static checks for kinetic SBML models. It finds rate laws that divide by zero,
reactions that can drive species negative, zero-order reactions that are
probably typos, and undefined or unused symbols. It reports each problem as a
finding, with a concrete witness where one exists.

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
