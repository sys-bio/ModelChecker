//! Every finding code, with a title and an explanation. Codes are stable:
//! never renumber or reuse one.
//!
//! Prefixes: `RL` rate laws, `SY` symbols and values, `NW` network,
//! `RU` rules and events, `UN` units, `MD` metadata.

pub struct Code {
    pub code: &'static str,
    pub title: &'static str,
    pub help: &'static str,
}

pub const RL_SINGULAR: &str = "RL001";
pub const RL_POSSIBLY_SINGULAR: &str = "RL002";
pub const RL_NEGATIVITY: &str = "RL003";
pub const RL_POSSIBLE_NEGATIVITY: &str = "RL004";
pub const RL_ZERO_ORDER: &str = "RL005";
pub const RL_MISSING_REACTANT: &str = "RL006";
pub const RL_CONSTANT_SOURCE: &str = "RL007";
pub const RL_NO_KINETIC_LAW: &str = "RL008";
pub const RL_UNDECLARED_SPECIES: &str = "RL009";
pub const RL_UNUSED_MODIFIER: &str = "RL010";
pub const RL_EMPTY_REACTION: &str = "RL011";
pub const RL_SINGULAR_INITIAL: &str = "RL012";
pub const RL_MASS_ACTION_ORDER: &str = "RL013";

pub const SY_UNDEFINED: &str = "SY001";
pub const SY_UNUSED_PARAMETER: &str = "SY002";
pub const SY_UNUSED_LOCAL: &str = "SY003";
pub const SY_UNUSED_FUNCTION: &str = "SY004";
pub const SY_SHADOWED: &str = "SY005";
pub const SY_NO_VALUE: &str = "SY006";
pub const SY_UNSUPPORTED_MATH: &str = "SY007";

pub const UN_RATE_UNITS: &str = "UN001";
pub const UN_DIMENSION_MISMATCH: &str = "UN002";
pub const UN_SCALE_MISMATCH: &str = "UN003";
pub const UN_FUNCTION_ARGUMENT: &str = "UN004";
pub const UN_MISSING_VOLUME: &str = "UN005";
pub const UN_VOLUME_HEURISTIC: &str = "UN006";
pub const UN_ASSIGNMENT: &str = "UN007";
pub const UN_COVERAGE: &str = "UN008";

pub const SO_UNKNOWN: &str = "SO001";
pub const SO_OBSOLETE: &str = "SO002";
pub const SO_WRONG_BRANCH: &str = "SO003";
pub const SO_IRREVERSIBLE_NEGATIVE: &str = "SO004";
pub const SO_REVERSIBLE_NEVER_NEGATIVE: &str = "SO005";
pub const SO_NOT_MASS_ACTION: &str = "SO006";
pub const SO_MASS_ACTION_MISMATCH: &str = "SO007";
pub const SO_SHAPE_MISMATCH: &str = "SO008";

pub const RU_RULE_LOOP: &str = "RU001";
pub const RU_INITIAL_LOOP: &str = "RU002";
pub const RU_SET_TWICE: &str = "RU003";
pub const RU_CONSTANT_CHANGED: &str = "RU004";
pub const RU_ALGEBRAIC: &str = "RU005";
pub const RU_TRIGGER_TRUE_AT_START: &str = "RU006";
pub const RU_EVENT_NEVER: &str = "RU007";
pub const RU_EVENT_EMPTY: &str = "RU008";
pub const RU_NEGATIVE_DELAY: &str = "RU009";

pub const NW_ONLY_PRODUCED: &str = "NW001";
pub const NW_ONLY_CONSUMED: &str = "NW002";
pub const NW_UNUSED_SPECIES: &str = "NW003";
pub const NW_DISCONNECTED: &str = "NW004";
pub const NW_CONSERVED: &str = "NW005";
pub const NW_NEVER_PRESENT: &str = "NW006";
pub const NW_NO_EFFECT: &str = "NW007";

pub const TD_HALDANE: &str = "TD001";
pub const TD_WEGSCHEIDER: &str = "TD002";

pub const AN_ERROR: &str = "AN001";
pub const AN_WARNING: &str = "AN002";

pub const SB_READ: &str = "SB001";
pub const SB_READ_WARNING: &str = "SB002";
pub const SB_STRUCTURE: &str = "SB003";
pub const SB_NO_MODEL: &str = "SB004";

pub const ALL: &[Code] = &[
    Code {
        code: RU_RULE_LOOP,
        title: "Rules in a loop",
        help: "Assignment rules depend on each other in a circle (in Level 3, possibly through \
               a reaction's rate), so none can be computed. SBML does not allow this.",
    },
    Code {
        code: RU_INITIAL_LOOP,
        title: "Initial values in a loop",
        help: "Initial assignments (with assignment rules) depend on each other in a circle, \
               so the initial values cannot be worked out. SBML does not allow this.",
    },
    Code {
        code: RU_SET_TWICE,
        title: "Set more than once",
        help: "A variable is set in two ways that conflict: two rules, an assignment rule and \
               an initial assignment, an assignment rule and an event, or reactions and a \
               rule on a species that is not a boundary species. SBML does not allow this.",
    },
    Code {
        code: RU_CONSTANT_CHANGED,
        title: "Constant changed",
        help: "Something marked constant is set by a rule or changed by an event.",
    },
    Code {
        code: RU_ALGEBRAIC,
        title: "Algebraic rules",
        help: "The model uses algebraic rules (constraints of the form 0 = f(...)). Many \
               simulators do not support them.",
    },
    Code {
        code: RU_TRIGGER_TRUE_AT_START,
        title: "Trigger true at the start",
        help: "An event's trigger is already true at the initial values. Depending on the \
               trigger's initialValue, the event either fires at once at time 0 or does not \
               fire until the trigger goes false and then true again; either may be a \
               surprise.",
    },
    Code {
        code: RU_EVENT_NEVER,
        title: "Event never fires",
        help: "The event's trigger depends on nothing that changes during a simulation, so it \
               can never become true.",
    },
    Code {
        code: RU_EVENT_EMPTY,
        title: "Event does nothing",
        help: "The event has no event assignments.",
    },
    Code {
        code: RU_NEGATIVE_DELAY,
        title: "Negative delay",
        help: "An event's delay is negative at the start; delays must not be negative.",
    },
    Code {
        code: NW_ONLY_PRODUCED,
        title: "Only produced",
        help: "No reaction consumes the species, so it accumulates. Often intended (a model \
               output); worth a look if other parts of the model use it.",
    },
    Code {
        code: NW_ONLY_CONSUMED,
        title: "Only consumed",
        help: "No reaction produces the species, so it runs down. Often intended (a substrate \
               pool).",
    },
    Code {
        code: NW_UNUSED_SPECIES,
        title: "Unused species",
        help: "The species takes part in no reaction, rule or event, and no expression uses it.",
    },
    Code {
        code: NW_DISCONNECTED,
        title: "Disconnected network",
        help: "The reactions fall into separate parts that never affect each other: no shared \
               species, and no rate law, rule or event linking them. Sometimes intended; \
               sometimes a missing link or a misspelt species.",
    },
    Code {
        code: NW_CONSERVED,
        title: "Conserved totals",
        help: "Conservation laws: weighted sums of species that no reaction changes, from the \
               left null space of the stoichiometry matrix. Useful to know; steady-state \
               solvers need them.",
    },
    Code {
        code: NW_NEVER_PRESENT,
        title: "Species that can never appear",
        help: "These species start at zero, and every reaction that would produce one of them \
               needs one of them as a reactant (an empty siphon), so they stay at zero and the \
               reactions consuming them never run. Usually a missing initial value.",
    },
    Code {
        code: NW_NO_EFFECT,
        title: "Reaction changes nothing",
        help: "The reaction changes no variable species: all its species are boundary or \
               constant, or appear equally on both sides. It only carries a flux.",
    },
    Code {
        code: TD_HALDANE,
        title: "Equilibrium not constant (Haldane)",
        help: "A reversible rate law is zero at equilibrium, which thermodynamics fixes at \
               one value of the mass-action ratio: the equilibrium constant. This rate law is \
               zero at different ratios depending on other variables, such as a modifier, so \
               it does not satisfy the Haldane relationship. Many published models are not \
               thermodynamically consistent, which may not matter for their purpose.",
    },
    Code {
        code: TD_WEGSCHEIDER,
        title: "Cycle of equilibria inconsistent (Wegscheider)",
        help: "Reversible reactions forming a closed cycle must have equilibrium constants \
               whose product around the cycle is 1; otherwise the cycle would turn forever at \
               equilibrium. Many published models are not thermodynamically consistent, which \
               may not matter for their purpose.",
    },
    Code {
        code: AN_ERROR,
        title: "Antimony error",
        help: "The edited Antimony text cannot be read. Fix it and re-check; the other \
               checks run once it reads.",
    },
    Code {
        code: AN_WARNING,
        title: "Antimony warning",
        help: "Something in the edited Antimony text is suspicious, such as a name that is \
               used but never declared, which is often a typo.",
    },
    Code {
        code: RL_MASS_ACTION_ORDER,
        title: "Mass-action order differs from stoichiometry",
        help: "The rate law has the form of mass action (a constant times powers of \
               species, minus a similar reverse term), but a species appears to a power \
               other than its stoichiometry, as in k*A for 2 A -> B. For an elementary \
               reaction this is a mistake; for a non-elementary one it may be intended.",
    },
    Code {
        code: SO_UNKNOWN,
        title: "Unknown SBO term",
        help: "The sboTerm is not a valid SBO id, or not a term in the Systems Biology \
               Ontology.",
    },
    Code {
        code: SO_OBSOLETE,
        title: "Obsolete SBO term",
        help: "The SBO term has been made obsolete in the ontology; a current term should be \
               used instead.",
    },
    Code {
        code: SO_WRONG_BRANCH,
        title: "SBO term of the wrong kind",
        help: "Each SBML element takes SBO terms from one branch of the ontology: rate laws \
               for kinetic laws, material entities for species and compartments, \
               participant roles for reactants and products, and so on. This term is from \
               another branch.",
    },
    Code {
        code: SO_IRREVERSIBLE_NEGATIVE,
        title: "Irreversible term, negative rate",
        help: "The kinetic law's SBO term describes an irreversible rate law, but the rate \
               can be negative, so the reaction can run backwards.",
    },
    Code {
        code: SO_REVERSIBLE_NEVER_NEGATIVE,
        title: "Reversible term, rate never negative",
        help: "The kinetic law's SBO term describes a reversible rate law, but the rate can \
               never be negative.",
    },
    Code {
        code: SO_NOT_MASS_ACTION,
        title: "Mass-action term, not mass action",
        help: "The kinetic law's SBO term is a mass-action rate law, but the rate law is not \
               a sum of terms that are each a constant times powers of species.",
    },
    Code {
        code: SO_MASS_ACTION_MISMATCH,
        title: "Mass-action term does not fit",
        help: "The kinetic law's SBO term gives orders and numbers of reactants and products \
               that differ from the rate law, or the rate law's terms do not match the \
               reaction's stoichiometry.",
    },
    Code {
        code: SO_SHAPE_MISMATCH,
        title: "Enzymatic or Hill term does not fit",
        help: "The kinetic law's SBO term names an enzymatic or Hill-type rate law, but the \
               rate law lacks that shape (saturation in the substrate, a species raised to a \
               power, a modifier) or the reaction has a different number of reactants or \
               products.",
    },
    Code {
        code: UN_RATE_UNITS,
        title: "Rate law units wrong",
        help: "A kinetic law must give extent per time (in Level 2, substance per time), but                the units of this rate law, worked out from its symbols' declared units, are                different.",
    },
    Code {
        code: UN_DIMENSION_MISMATCH,
        title: "Units do not match",
        help: "An expression adds, subtracts or compares quantities of different dimensions,                for example a concentration and a rate constant.",
    },
    Code {
        code: UN_SCALE_MISMATCH,
        title: "Units differ in scale",
        help: "An expression adds, subtracts or compares quantities of the same kind in                different units, for example mM and M, without converting.",
    },
    Code {
        code: UN_FUNCTION_ARGUMENT,
        title: "Argument not dimensionless",
        help: "exp, ln, log, trigonometric functions and exponents need dimensionless                arguments; this one has units.",
    },
    Code {
        code: UN_MISSING_VOLUME,
        title: "Missing volume factor",
        help: "The rate law is in concentration per time, short of amount per time by exactly                a compartment's size. Multiply it by the compartment.",
    },
    Code {
        code: UN_VOLUME_HEURISTIC,
        title: "Possibly missing volume factor",
        help: "The model's units are not enough to tell, but the rate law uses concentrations                of species in a compartment whose size is not 1 and never mentions a                compartment. Kinetic laws give amount per time, so this is usually missing a                factor of the compartment size.",
    },
    Code {
        code: UN_ASSIGNMENT,
        title: "Rule units wrong",
        help: "An assignment rule, rate rule, initial assignment or event assignment gives a                value whose units differ from its variable's (per time, for a rate rule).",
    },
    Code {
        code: UN_COVERAGE,
        title: "Units not fully checked",
        help: "How much of the model's units could not be checked, and why: usually                parameters with no declared units.",
    },
    Code {
        code: SB_READ,
        title: "Cannot read SBML",
        help: "The file is not well-formed XML, or not SBML Level 2 or 3.",
    },
    Code {
        code: SB_READ_WARNING,
        title: "Reader warning",
        help: "Something unusual was found while reading the file.",
    },
    Code {
        code: SB_STRUCTURE,
        title: "Structural problem",
        help: "Duplicate ids, or a reference to an element that does not exist.",
    },
    Code {
        code: SB_NO_MODEL,
        title: "No model",
        help: "The document contains no model.",
    },
    Code {
        code: RL_SINGULAR,
        title: "Rate cannot be calculated at some values",
        help: "For some values the model could reach, the rate law cannot be calculated: it \
               would divide by zero, take the logarithm of zero or a negative number, take the \
               square root of a negative number, or similar. The witness gives such values. \
               Whether a simulation actually gets there depends on how the model behaves; \
               some models (power-law models, for example) are written this way on purpose.",
    },
    Code {
        code: RL_SINGULAR_INITIAL,
        title: "Rate cannot be calculated at the start",
        help: "At the model's starting values, the rate law cannot be calculated: it divides \
               by zero, takes the logarithm of zero or a negative number, or similar. A \
               simulation will fail immediately.",
    },
    Code {
        code: RL_POSSIBLY_SINGULAR,
        title: "Rate might not be calculable",
        help: "The checker could not rule out that the rate law divides by zero (or similar) \
               for some values of the species, but found no actual values where it does. \
               Often a false alarm for complicated rate laws; worth a look.",
    },
    Code {
        code: RL_NEGATIVITY,
        title: "Species can go negative",
        help: "With the species at zero, the reaction still consumes it, so a simulation can \
               drive it negative. For a net reactant the rate must be <= 0 when the species \
               is zero; for a net product it must be >= 0.",
    },
    Code {
        code: RL_POSSIBLE_NEGATIVITY,
        title: "Species may go negative",
        help: "Could not show that the reaction stops consuming the species when it reaches \
               zero, but no concrete counterexample was found.",
    },
    Code {
        code: RL_ZERO_ORDER,
        title: "Zero-order reaction",
        help: "The reaction consumes species but its rate depends on none of them. This can \
               be intended (saturated enzyme), but it is often a typo, and the reactants can \
               go negative.",
    },
    Code {
        code: RL_MISSING_REACTANT,
        title: "Reactant missing from rate law",
        help: "The reaction consumes this species but its rate does not depend on it, so the \
               species can go negative.",
    },
    Code {
        code: RL_CONSTANT_SOURCE,
        title: "Constant-rate source",
        help: "The reaction consumes no variable species and its rate depends on no species. \
               Usually an intended inflow.",
    },
    Code {
        code: RL_NO_KINETIC_LAW,
        title: "No kinetic law",
        help: "The reaction has no rate law, so it cannot be simulated.",
    },
    Code {
        code: RL_UNDECLARED_SPECIES,
        title: "Species in rate law not declared",
        help: "The rate law refers to a species that is not a reactant, product or modifier \
               of the reaction. Declare it as a modifier.",
    },
    Code {
        code: RL_UNUSED_MODIFIER,
        title: "Modifier not used",
        help: "The species is declared as a modifier but the rate law does not depend on it.",
    },
    Code {
        code: RL_EMPTY_REACTION,
        title: "Empty reaction",
        help: "The reaction has no reactants and no products.",
    },
    Code {
        code: SY_UNDEFINED,
        title: "Undefined symbol",
        help: "An expression refers to an id that is not defined anywhere in the model.",
    },
    Code {
        code: SY_UNUSED_PARAMETER,
        title: "Unused parameter",
        help: "The parameter is defined but no expression refers to it and nothing sets it.",
    },
    Code {
        code: SY_UNUSED_LOCAL,
        title: "Unused local parameter",
        help: "The local parameter is not used in its reaction's rate law.",
    },
    Code {
        code: SY_UNUSED_FUNCTION,
        title: "Unused function",
        help: "The function definition is never called.",
    },
    Code {
        code: SY_SHADOWED,
        title: "Local parameter shadows global",
        help: "A local parameter has the same id as a global symbol, which it hides inside \
               the rate law. Legal, but easy to misread.",
    },
    Code {
        code: SY_NO_VALUE,
        title: "No value",
        help: "The symbol has no value, no initial assignment and no rule that sets it.",
    },
    Code {
        code: SY_UNSUPPORTED_MATH,
        title: "Math not analysed",
        help: "The expression contains MathML the checker does not understand, so it was \
               skipped.",
    },
];

pub fn lookup(code: &str) -> Option<&'static Code> {
    ALL.iter().find(|c| c.code == code)
}
