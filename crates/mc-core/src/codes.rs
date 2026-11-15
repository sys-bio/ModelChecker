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

pub const SB_READ: &str = "SB001";
pub const SB_READ_WARNING: &str = "SB002";
pub const SB_STRUCTURE: &str = "SB003";
pub const SB_NO_MODEL: &str = "SB004";

pub const ALL: &[Code] = &[
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
        title: "Singular rate law",
        help: "A division by zero, log of a non-positive number, root of a negative number, \
               or similar, occurs at a state with all species non-negative (and no conserved \
               total forced to zero). Whether a simulation gets there depends on the \
               dynamics; power-law models are singular at zero by design. The witness gives \
               a point where it happens.",
    },
    Code {
        code: RL_SINGULAR_INITIAL,
        title: "Singular at initial values",
        help: "A division by zero, log of a non-positive number, or similar, occurs at the \
               model's initial values (or for all values), so a simulation fails at once.",
    },
    Code {
        code: RL_POSSIBLY_SINGULAR,
        title: "Possibly singular rate law",
        help: "Interval analysis could not rule out a division by zero (or similar) for \
               non-negative species values, but no concrete point was found. Often a false \
               alarm for complex expressions; worth a look.",
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
