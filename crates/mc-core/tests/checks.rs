//! Each check on small models written in Antimony.

use mc_core::codes::*;
use mc_core::{Report, Severity};

fn doc(ant: &str) -> sbml_io::Document {
    let loaded = antimony::parse(ant);
    assert!(!loaded.has_errors(), "{:?}", loaded.diagnostics);
    let flat = antimony::flatten::flatten(&loaded.document);
    antimony::sbml::to_sbml(&flat)
}

fn run(ant: &str) -> Report {
    mc_core::check(&sbml_io::write(&doc(ant)))
}

fn codes(r: &Report) -> Vec<&'static str> {
    r.findings.iter().map(|f| f.code).collect()
}

fn has(r: &Report, code: &str) -> bool {
    r.findings.iter().any(|f| f.code == code)
}

/// Fails with the whole report, so a failing test shows what was found.
macro_rules! assert_has {
    ($r:expr, $code:expr) => {
        assert!(has(&$r, $code), "expected {} in:\n{}", $code, show(&$r))
    };
}
macro_rules! assert_not {
    ($r:expr, $code:expr) => {
        assert!(!has(&$r, $code), "did not expect {} in:\n{}", $code, show(&$r))
    };
}

fn show(r: &Report) -> String {
    r.findings.iter().map(|f| format!("  {f}\n")).collect()
}

fn rate_law_codes(r: &Report) -> Vec<&'static str> {
    codes(r).into_iter().filter(|c| c.starts_with("RL")).collect()
}

#[test]
fn mass_action_is_clean() {
    let r = run("J1: S1 => S2; k1*S1\nS1 = 10; S2 = 0; k1 = 0.1");
    assert!(rate_law_codes(&r).is_empty(), "{}", show(&r));
}

#[test]
fn reversible_mass_action_is_clean() {
    let r = run("J1: A -> B; k1*A - k2*B\nA = 1; B = 0; k1 = 1; k2 = 0.5");
    assert!(rate_law_codes(&r).is_empty(), "{}", show(&r));
}

#[test]
fn reversible_michaelis_menten_is_clean() {
    let r = run(
        "J1: S -> P; (Vf*S/Ks - Vr*P/Kp)/(1 + S/Ks + P/Kp)\n\
         S = 1; P = 0; Vf = 1; Vr = 0.5; Ks = 0.2; Kp = 0.3",
    );
    assert!(rate_law_codes(&r).is_empty(), "{}", show(&r));
}

#[test]
fn zero_order() {
    let r = run("J1: S1 => S2; k1\nS1 = 10; S2 = 0; k1 = 0.1");
    assert_has!(r, RL_ZERO_ORDER);
    // Reported once, not again as a positivity failure.
    assert_not!(r, RL_NEGATIVITY);
}

#[test]
fn missing_reactant() {
    let r = run("J1: A + B => C; k*A\nA = 1; B = 1; C = 0; k = 1");
    assert_has!(r, RL_MISSING_REACTANT);
    assert!(r.findings.iter().any(|f| f.message.starts_with("The reaction consumes B")));
    assert_not!(r, RL_ZERO_ORDER);
}

#[test]
fn constant_source() {
    let r = run("J1: => S1; v0\nS1 = 0; v0 = 1");
    assert_has!(r, RL_CONSTANT_SOURCE);
    assert!(r.findings.iter().all(|f| f.severity != Severity::Error), "{}", show(&r));
}

#[test]
fn boundary_source() {
    let r = run("J1: $X0 => S1; v0\nX0 = 1; S1 = 0; v0 = 1");
    assert_has!(r, RL_CONSTANT_SOURCE);
    assert_not!(r, RL_ZERO_ORDER);
}

#[test]
fn reactant_consumed_at_zero() {
    let r = run("J1: S => P; k*(S + 1)\nS = 1; P = 0; k = 1");
    assert_has!(r, RL_NEGATIVITY);
    let f = r.findings.iter().find(|f| f.code == RL_NEGATIVITY).unwrap();
    assert!(f.witness.as_deref().unwrap().starts_with("S = 0"), "{f}");
}

#[test]
fn product_consumed_at_zero() {
    // At P = 0 and small S the rate is negative: the reaction runs
    // backwards and consumes P, which is not there.
    let r = run("J1: S -> P; k1*S - k2\nS = 1; P = 0; k1 = 1; k2 = 0.5");
    let f = r
        .findings
        .iter()
        .find(|f| f.code == RL_NEGATIVITY)
        .unwrap_or_else(|| panic!("{}", show(&r)));
    assert!(f.witness.as_deref().unwrap().starts_with("P = 0"), "{f}");
}

#[test]
fn division_by_difference() {
    let r = run("J1: S -> P; k*S/(S - P)\nS = 1; P = 0; k = 1");
    assert_has!(r, RL_SINGULAR);
}

#[test]
fn zero_michaelis_constant() {
    let r = run("J1: S => P; Vm*S/(Km + S)\nS = 1; P = 0; Vm = 1; Km = 0");
    assert_has!(r, RL_SINGULAR);
}

#[test]
fn michaelis_menten_is_clean() {
    let r = run("J1: S => P; Vm*S/(Km + S)\nS = 1; P = 0; Vm = 1; Km = 0.5");
    assert!(rate_law_codes(&r).is_empty(), "{}", show(&r));
}

#[test]
fn hill_is_clean() {
    let r = run("J1: S => P; Vm*S^h/(K^h + S^h)\nS = 1; P = 0; Vm = 1; K = 0.5; h = 2.5");
    assert!(rate_law_codes(&r).is_empty(), "{}", show(&r));
}

#[test]
fn guarded_division_is_not_reported() {
    let r = run("J1: S => P; piecewise(0, S == 0, k*S*S/S)\nS = 1; P = 0; k = 1");
    assert_not!(r, RL_SINGULAR);
    assert_not!(r, RL_POSSIBLY_SINGULAR);
}

#[test]
fn log_of_species() {
    let r = run("J1: S => P; k*S*ln(S)\nS = 1; P = 0; k = 1");
    assert_has!(r, RL_SINGULAR);
}

#[test]
fn through_function_definition() {
    let clean = run(
        "function MM(S, Vm, Km)\n  Vm*S/(Km + S)\nend\n\
         J1: S => P; MM(S, Vm, Km)\nS = 1; P = 0; Vm = 1; Km = 0.5",
    );
    assert!(rate_law_codes(&clean).is_empty(), "{}", show(&clean));
    let zero = run(
        "function f(k)\n  2*k\nend\n\
         J1: S => P; f(k)\nS = 1; P = 0; k = 1",
    );
    assert_has!(zero, RL_ZERO_ORDER);
}

#[test]
fn through_assignment_rule() {
    let r = run("x := 2*S\nJ1: S => P; k*x\nS = 1; P = 0; k = 1");
    assert_not!(r, RL_ZERO_ORDER);
    assert_not!(r, RL_NEGATIVITY);
}

#[test]
fn singular_assignment_rule_reported_once() {
    let r = run(
        "x := 1/(S - P)\nJ1: S -> P; k*S*x\nJ2: S -> P; k*S*x\n\
         S = 1; P = 0; k = 1",
    );
    let n = r.findings.iter().filter(|f| f.code == RL_SINGULAR).count();
    assert_eq!(n, 1, "{}", show(&r));
}

#[test]
fn unused_parameter() {
    let r = run("J1: S => P; k*S\nS = 1; P = 0; k = 1; k2 = 3");
    assert_has!(r, SY_UNUSED_PARAMETER);
}

#[test]
fn undeclared_species_and_undefined_symbol() {
    let mut d = doc("species E\nJ1: S => P; k*S*E\nS = 1; P = 0; E = 1; k = 1");
    let m = d.model.as_mut().unwrap();
    m.reactions[0].modifiers.clear();
    let r = mc_core::check(&sbml_io::write(&d));
    assert_has!(r, RL_UNDECLARED_SPECIES);

    let mut d = doc("J1: S => P; k*S\nS = 1; P = 0; k = 1");
    let m = d.model.as_mut().unwrap();
    m.reactions[0].kinetic_law.as_mut().unwrap().math = Some(sbml_io::Math::Apply(
        sbml_io::Op::Times,
        vec![sbml_io::Math::ident("zz"), sbml_io::Math::ident("S")],
    ));
    let r = mc_core::check(&sbml_io::write(&d));
    assert_has!(r, SY_UNDEFINED);
}

#[test]
fn unreadable_file() {
    let r = mc_core::check("not xml");
    assert_has!(r, SB_READ);
}

#[test]
fn conserved_total_in_denominator_is_not_singular() {
    // E + ES is conserved and starts positive, so E + ES never vanishes.
    let r = run(
        "J1: E + S -> ES; k1*E*S - k2*ES\n\
         J2: ES => E + P; k3*ES/(E + ES)\n\
         E = 1; S = 10; ES = 0; P = 0; k1 = 1; k2 = 1; k3 = 1",
    );
    assert_not!(r, RL_SINGULAR);
    assert_not!(r, RL_SINGULAR_INITIAL);
}

#[test]
fn singular_at_initial_values_is_an_error() {
    let r = run("J1: S => P; k*S/P\nS = 1; P = 0; k = 1");
    let f = r
        .findings
        .iter()
        .find(|f| f.code == RL_SINGULAR_INITIAL)
        .unwrap_or_else(|| panic!("{}", show(&r)));
    assert_eq!(f.severity, Severity::Error);
}

#[test]
fn negative_species_is_not_held_to_positivity() {
    // V is a voltage: it starts negative, so it is not a concentration.
    let r = run("J1: => V; gL*(EL - V)\nV = -65; gL = 0.1; EL = -70");
    assert_not!(r, RL_NEGATIVITY);
    assert_not!(r, RL_POSSIBLE_NEGATIVITY);
}

// ---------------------------------------------------------------- units

fn unit_def(id: &str, parts: &[(&str, f64, i32)]) -> sbml_io::UnitDefinition {
    sbml_io::UnitDefinition {
        id: id.into(),
        units: parts
            .iter()
            .map(|(kind, exponent, scale)| sbml_io::Unit {
                kind: kind.to_string(),
                exponent: Some(*exponent),
                scale: Some(*scale),
                multiplier: Some(1.0),
                ..Default::default()
            })
            .collect(),
        ..Default::default()
    }
}

/// `ant` with mole, litre and second as the model's units, and units for
/// the named parameters: `per_s`, `M` (mol/l), `mM`, `M_per_s`.
fn run_units(ant: &str, params: &[(&str, &str)]) -> Report {
    let mut d = doc(ant);
    let m = d.model.as_mut().unwrap();
    m.substance_units = Some("mole".into());
    m.extent_units = Some("mole".into());
    m.time_units = Some("second".into());
    m.volume_units = Some("litre".into());
    m.unit_definitions = vec![
        unit_def("per_s", &[("second", -1.0, 0)]),
        unit_def("M", &[("mole", 1.0, 0), ("litre", -1.0, 0)]),
        unit_def("mM", &[("mole", 1.0, -3), ("litre", -1.0, 0)]),
        unit_def("M_per_s", &[("mole", 1.0, 0), ("litre", -1.0, 0), ("second", -1.0, 0)]),
    ];
    for (p, u) in params {
        let par = m.parameters.iter_mut().find(|x| x.id == *p).unwrap();
        par.units = Some(u.to_string());
    }
    mc_core::check(&sbml_io::write(&d))
}

fn unit_codes(r: &Report) -> Vec<&'static str> {
    codes(r).into_iter().filter(|c| c.starts_with("UN")).collect()
}

#[test]
fn units_consistent_michaelis_menten() {
    let r = run_units(
        "compartment cell = 2\nspecies S in cell, P in cell\n\
         J1: S => P; cell*Vm*S/(Km + S)\nS = 1; P = 0; Vm = 1; Km = 0.5",
        &[("Vm", "M_per_s"), ("Km", "M")],
    );
    assert!(unit_codes(&r).is_empty(), "{}", show(&r));
}

#[test]
fn units_missing_volume_factor() {
    let r = run_units(
        "compartment cell = 2\nspecies S in cell, P in cell\n\
         J1: S => P; Vm*S/(Km + S)\nS = 1; P = 0; Vm = 1; Km = 0.5",
        &[("Vm", "M_per_s"), ("Km", "M")],
    );
    assert_has!(r, UN_MISSING_VOLUME);
    assert_not!(r, UN_RATE_UNITS);
}

#[test]
fn units_wrong_rate() {
    // k*S in a 1-litre compartment with k per second gives mol/(l*s), short
    // of the volume; k in M_per_s gives M^2/s: just wrong.
    let r = run_units(
        "compartment cell = 1\nspecies S in cell, P in cell\n\
         J1: S => P; cell*k*S\nS = 1; P = 0; k = 1",
        &[("k", "M_per_s")],
    );
    assert_has!(r, UN_RATE_UNITS);
}

#[test]
fn units_adding_different_dimensions() {
    let r = run_units(
        "compartment cell = 1\nspecies S in cell, P in cell\n\
         J1: S => P; cell*k*(S + k)\nS = 1; P = 0; k = 1",
        &[("k", "per_s")],
    );
    assert_has!(r, UN_DIMENSION_MISMATCH);
}

#[test]
fn units_adding_different_scales() {
    let r = run_units(
        "compartment cell = 1\nspecies S in cell, P in cell\n\
         J1: S => P; cell*Vm*S/(Km + S)\nS = 1; P = 0; Vm = 1; Km = 0.5",
        &[("Vm", "M_per_s"), ("Km", "mM")],
    );
    assert_has!(r, UN_SCALE_MISMATCH);
}

#[test]
fn units_exp_of_concentration() {
    let r = run_units(
        "compartment cell = 1\nspecies S in cell, P in cell\n\
         J1: S => P; cell*k*S*exp(S)\nS = 1; P = 0; k = 1",
        &[("k", "per_s")],
    );
    assert_has!(r, UN_FUNCTION_ARGUMENT);
}

#[test]
fn units_numbers_adapt() {
    // 1 + S/Ki: the 1 takes the units of S/Ki (dimensionless).
    let r = run_units(
        "compartment cell = 1\nspecies S in cell, P in cell\n\
         J1: S => P; cell*k*S/(1 + S/Ki)\nS = 1; P = 0; k = 1; Ki = 2",
        &[("k", "per_s"), ("Ki", "M")],
    );
    assert!(unit_codes(&r).is_empty(), "{}", show(&r));
}

#[test]
fn units_assignment_rule() {
    let r = run_units(
        "compartment cell = 1\nspecies S in cell\nx := S*k\nS = 1; k = 1\nx = 0",
        &[("k", "per_s"), ("x", "M")],
    );
    assert_has!(r, UN_ASSIGNMENT);
}

#[test]
fn units_coverage_without_units() {
    let r = run("J1: S => P; k*S\nS = 1; P = 0; k = 1");
    assert_has!(r, UN_COVERAGE);
    assert_not!(r, UN_RATE_UNITS);
}

#[test]
fn volume_heuristic_without_units() {
    let r = run("compartment cell = 2\nspecies S in cell, P in cell\nJ1: S => P; k*S\nS = 1; P = 0; k = 1");
    assert_has!(r, UN_VOLUME_HEURISTIC);
    let r = run("compartment cell = 2\nspecies S in cell, P in cell\nJ1: S => P; cell*k*S\nS = 1; P = 0; k = 1");
    assert_not!(r, UN_VOLUME_HEURISTIC);
    let r = run("compartment cell = 1\nspecies S in cell, P in cell\nJ1: S => P; k*S\nS = 1; P = 0; k = 1");
    assert_not!(r, UN_VOLUME_HEURISTIC);
}

#[test]
fn repeated_unit_findings_are_merged() {
    // Three rate laws missing the same volume factor: one finding.
    let r = run_units(
        "compartment cell = 2\nspecies A in cell, B in cell, C in cell, D in cell\n\
         J1: A => B; k*A\nJ2: B => C; k*B\nJ3: C => D; k*C\n\
         A = 1; B = 0; C = 0; D = 0; k = 1",
        &[("k", "per_s")],
    );
    let v: Vec<_> = r.findings.iter().filter(|f| f.code == UN_MISSING_VOLUME).collect();
    assert_eq!(v.len(), 1, "{}", show(&r));
    assert_eq!(v[0].element, mc_core::Element::Model);
    assert!(v[0].message.starts_with("3 rate laws"), "{}", v[0].message);
    assert!(v[0].message.contains("J1, J2, J3"), "{}", v[0].message);

    // Two are reported one by one.
    let r = run_units(
        "compartment cell = 2\nspecies A in cell, B in cell, C in cell\n\
         J1: A => B; k*A\nJ2: B => C; k*B\nA = 1; B = 0; C = 0; k = 1",
        &[("k", "per_s")],
    );
    let n = r.findings.iter().filter(|f| f.code == UN_MISSING_VOLUME).count();
    assert_eq!(n, 2, "{}", show(&r));
}

// ---------------------------------------------------------------- SBO

/// `ant` with the first reaction's kinetic law given `sbo`.
fn run_sbo(ant: &str, sbo: &str) -> Report {
    let mut d = doc(ant);
    let m = d.model.as_mut().unwrap();
    m.reactions[0].kinetic_law.as_mut().unwrap().meta.sbo_term = Some(sbo.into());
    mc_core::check(&sbml_io::write(&d))
}

fn sbo_codes(r: &Report) -> Vec<&'static str> {
    codes(r).into_iter().filter(|c| c.starts_with("SO")).collect()
}

const MASS_ACTION_1: &str = "J1: S => P; k*S\nS = 1; P = 0; k = 1";

#[test]
fn sbo_first_order_mass_action_fits() {
    let r = run_sbo(MASS_ACTION_1, "SBO:0000049");
    assert!(sbo_codes(&r).is_empty(), "{}", show(&r));
}

#[test]
fn sbo_mass_action_order_wrong() {
    let r = run_sbo("species E\nJ1: S => P; k*S*E\nS = 1; P = 0; E = 1; k = 1", "SBO:0000049");
    assert_has!(r, SO_MASS_ACTION_MISMATCH);
}

#[test]
fn sbo_mass_action_against_stoichiometry() {
    // Second order, two reactants, but the reaction is 2 A -> C.
    let r = run_sbo("J1: 2 A => C; k*A*B\nA = 1; B = 1; C = 0; k = 1", "SBO:0000054");
    assert_has!(r, SO_MASS_ACTION_MISMATCH);
    let r = run_sbo("J1: A + B => C; k*A*B\nA = 1; B = 1; C = 0; k = 1", "SBO:0000054");
    assert!(sbo_codes(&r).is_empty(), "{}", show(&r));
}

#[test]
fn sbo_irreversible_but_negative() {
    let r = run_sbo("J1: S -> P; k1*S - k2*P\nS = 1; P = 0; k1 = 1; k2 = 1", "SBO:0000049");
    assert_has!(r, SO_IRREVERSIBLE_NEGATIVE);
}

#[test]
fn sbo_reversible_mass_action() {
    let r = run_sbo("J1: A -> B; k1*A - k2*B\nA = 1; B = 0; k1 = 1; k2 = 1", "SBO:0000080");
    assert!(sbo_codes(&r).is_empty(), "{}", show(&r));
    let r = run_sbo("J1: A -> B; k1*A\nA = 1; B = 0; k1 = 1", "SBO:0000080");
    assert_has!(r, SO_REVERSIBLE_NEVER_NEGATIVE);
    assert_has!(r, SO_MASS_ACTION_MISMATCH);
}

#[test]
fn sbo_mass_action_term_on_michaelis_menten() {
    let r = run_sbo("J1: S => P; Vm*S/(Km + S)\nS = 1; P = 0; Vm = 1; Km = 1", "SBO:0000049");
    assert_has!(r, SO_NOT_MASS_ACTION);
}

#[test]
fn sbo_michaelis_menten() {
    let mm = "J1: S => P; Vm*S/(Km + S)\nS = 1; P = 0; Vm = 1; Km = 1";
    let r = run_sbo(mm, "SBO:0000028");
    assert!(sbo_codes(&r).is_empty(), "{}", show(&r));
    let r = run_sbo(MASS_ACTION_1, "SBO:0000028");
    assert_has!(r, SO_SHAPE_MISMATCH);
    // Hill needs a power.
    let r = run_sbo(mm, "SBO:0000192");
    assert_has!(r, SO_SHAPE_MISMATCH);
    let r = run_sbo("J1: S => P; Vm*S^h/(K^h + S^h)\nS = 1; P = 0; Vm = 1; K = 1; h = 2", "SBO:0000192");
    assert!(sbo_codes(&r).is_empty(), "{}", show(&r));
}

#[test]
fn sbo_wrong_branch_and_unknown() {
    let r = run_sbo(MASS_ACTION_1, "SBO:0000176"); // biochemical reaction
    assert_has!(r, SO_WRONG_BRANCH);
    let r = run_sbo(MASS_ACTION_1, "SBO:9999999");
    assert_has!(r, SO_UNKNOWN);
    let mut d = doc(MASS_ACTION_1);
    d.model.as_mut().unwrap().species[0].meta.sbo_term = Some("SBO:0000049".into());
    let r = mc_core::check(&sbml_io::write(&d));
    assert_has!(r, SO_WRONG_BRANCH);
}

#[test]
fn sbo_mass_action_with_template_modifier() {
    // Translation from a template: first order in the modifier X.
    let r = run_sbo("species X\nJ1: => P; k*X\nX = 1; P = 0; k = 1", "SBO:0000049");
    assert_not!(r, SO_MASS_ACTION_MISMATCH);
}

// ------------------------------------------------- mass-action order

#[test]
fn mass_action_order_matches_stoichiometry() {
    let r = run("J1: 2 A => B; k*A^2\nA = 1; B = 0; k = 1");
    assert_not!(r, RL_MASS_ACTION_ORDER);
    let r = run("J1: A + B -> 2 C; k1*A*B - k2*C^2\nA = 1; B = 1; C = 0; k1 = 1; k2 = 1");
    assert_not!(r, RL_MASS_ACTION_ORDER);
    // Autocatalysis: X on both sides, first order in X.
    let r = run("J1: A + X => 2 X; k*A*X\nA = 1; X = 1; k = 1");
    assert_not!(r, RL_MASS_ACTION_ORDER);
}

#[test]
fn mass_action_order_differs() {
    let r = run("J1: 2 A => B; k*A\nA = 1; B = 0; k = 1");
    assert_has!(r, RL_MASS_ACTION_ORDER);
    let r = run("J1: A -> 2 B; k1*A - k2*B\nA = 1; B = 0; k1 = 1; k2 = 1");
    assert_has!(r, RL_MASS_ACTION_ORDER);
}

#[test]
fn mass_action_order_ignores_other_shapes() {
    // Michaelis-Menten is not mass action; a modifier-driven synthesis has
    // no reactant to compare.
    let r = run("J1: 2 A => B; Vm*A/(Km + A)\nA = 1; B = 0; Vm = 1; Km = 1");
    assert_not!(r, RL_MASS_ACTION_ORDER);
    let r = run("species E\nJ1: => P; k*E\nE = 1; P = 0; k = 1");
    assert_not!(r, RL_MASS_ACTION_ORDER);
}

#[test]
fn mass_action_order_skips_power_law() {
    let r = run("J1: A => B; a*A^g\nA = 1; B = 0; a = 1; g = 0.5");
    assert_not!(r, RL_MASS_ACTION_ORDER);
}

// ------------------------------------------------- thermodynamics

#[test]
fn haldane_consistent_laws() {
    let r = run(
        "J1: S -> P; (Vf*S/Ks - Vr*P/Kp)/(1 + S/Ks + P/Kp)\n\
         S = 1; P = 0; Vf = 1; Vr = 0.5; Ks = 0.2; Kp = 0.3",
    );
    assert_not!(r, TD_HALDANE);
    // An enzyme scaling both directions does not shift the equilibrium.
    let r = run("species E\nJ1: S -> P; E*(k1*S - k2*P)\nS = 1; P = 0; E = 1; k1 = 1; k2 = 2");
    assert_not!(r, TD_HALDANE);
}

#[test]
fn haldane_violated_by_a_modifier() {
    // E speeds up only the forward direction: it shifts the equilibrium.
    let r = run("species E\nJ1: S -> P; E*k1*S - k2*P\nS = 1; P = 0; E = 1; k1 = 1; k2 = 2");
    let f = r
        .findings
        .iter()
        .find(|f| f.code == TD_HALDANE)
        .unwrap_or_else(|| panic!("{}", show(&r)));
    assert!(f.message.contains("moves with E"), "{}", f.message);
}

#[test]
fn wegscheider_cycle() {
    let cycle = |k3r: &str| {
        run(&format!(
            "J1: A -> B; k1*A - k1r*B\nJ2: B -> C; k2*B - k2r*C\nJ3: C -> A; k3*C - k3r*A\n\
             A = 1; B = 1; C = 1; k1 = 2; k1r = 1; k2 = 3; k2r = 1; k3 = 1; k3r = {k3r}"
        ))
    };
    // K1*K2*K3 = 2 * 3 * (1/6) = 1: consistent.
    assert_not!(cycle("6"), TD_WEGSCHEIDER);
    // 2 * 3 * 1 = 6: the cycle would turn at equilibrium.
    let r = cycle("1");
    assert_has!(r, TD_WEGSCHEIDER);
    assert!(r.findings.iter().any(|f| f.message.contains("multiply to 6")), "{}", show(&r));
}

#[test]
fn wegscheider_ignores_paths_between_boundaries() {
    let r = run(
        "J1: $X -> A; k1*X - k1r*A\nJ2: A -> $Y; k2*A - k2r*Y\n\
         X = 1; Y = 1; A = 1; k1 = 2; k1r = 1; k2 = 3; k2r = 1",
    );
    assert_not!(r, TD_WEGSCHEIDER);
}

// ------------------------------------------------- network

#[test]
fn only_produced_and_only_consumed() {
    let r = run("J1: S => P; k*S\nS = 10; P = 0; k = 1");
    assert!(r.findings.iter().any(|f| f.code == NW_ONLY_PRODUCED && f.element == mc_core::Element::Species("P".into())), "{}", show(&r));
    assert!(r.findings.iter().any(|f| f.code == NW_ONLY_CONSUMED && f.element == mc_core::Element::Species("S".into())), "{}", show(&r));
}

#[test]
fn conserved_total() {
    let r = run(
        "J1: E + S -> ES; k1*E*S - k2*ES\nJ2: ES => E + P; k3*ES\n\
         E = 1; S = 10; ES = 0; P = 0; k1 = 1; k2 = 1; k3 = 1",
    );
    let f = r.findings.iter().find(|f| f.code == NW_CONSERVED).unwrap_or_else(|| panic!("{}", show(&r)));
    // E + ES is conserved with total 1; S + ES + P with total 10.
    assert!(f.message.contains("= 1\n") || f.message.ends_with("= 1"), "{}", f.message);
    assert!(f.message.contains("= 10"), "{}", f.message);
}

#[test]
fn disconnected_network() {
    let r = run("J1: A => B; k*A\nJ2: C => D; k*C\nA = 1; B = 0; C = 1; D = 0; k = 1");
    assert_has!(r, NW_DISCONNECTED);
    let r = run("J1: A => B; k*A\nJ2: B => C; k*B\nA = 1; B = 0; C = 0; k = 1");
    assert_not!(r, NW_DISCONNECTED);
}

#[test]
fn species_that_can_never_appear() {
    // X is autocatalytic and starts at zero: it can never appear.
    let r = run("J1: A + X => 2 X; k*A*X\nJ2: X => ; d*X\nA = 10; X = 0; k = 1; d = 0.1");
    let f = r.findings.iter().find(|f| f.code == NW_NEVER_PRESENT).unwrap_or_else(|| panic!("{}", show(&r)));
    assert!(f.message.starts_with("X starts at zero"), "{}", f.message);
    // With a seed of X it can.
    let r = run("J1: A + X => 2 X; k*A*X\nJ2: X => ; d*X\nA = 10; X = 0.01; k = 1; d = 0.1");
    assert_not!(r, NW_NEVER_PRESENT);
    // A chain from a source can always fill.
    let r = run("J0: => A; v\nJ1: A => B; k*A\nA = 0; B = 0; v = 1; k = 1");
    assert_not!(r, NW_NEVER_PRESENT);
}

#[test]
fn unused_species_and_no_effect_reaction() {
    let r = run("species Z\nJ1: S => P; k*S\nS = 1; P = 0; Z = 1; k = 1");
    assert_has!(r, NW_UNUSED_SPECIES);
    let r = run("J1: $X => $Y; k*X\nX = 1; Y = 0; k = 1");
    assert_has!(r, NW_NO_EFFECT);
}

// ------------------------------------------------- rules and events

#[test]
fn rules_in_a_loop() {
    // Antimony refuses circular rules, so make the loop in the SBML: y := 2*x.
    let mut d = doc("x := y + 1\ny := 2*S\nJ1: S => P; k*S*x\nS = 1; P = 0; k = 1");
    for rule in &mut d.model.as_mut().unwrap().rules {
        if let sbml_io::Rule::Assignment { variable, math, .. } = rule
            && variable == "y"
        {
            *math = Some(sbml_io::Math::Apply(
                sbml_io::Op::Times,
                vec![sbml_io::Math::int(2), sbml_io::Math::ident("x")],
            ));
        }
    }
    let r = mc_core::check(&sbml_io::write(&d));
    let f = r.findings.iter().find(|f| f.code == RU_RULE_LOOP).unwrap_or_else(|| panic!("{}", show(&r)));
    assert_eq!(f.severity, Severity::Error);
    assert!(f.message.contains("x depends on y") || f.message.contains("y depends on x"), "{}", f.message);
    let r = run("x := y + 1\ny := 2*S\nJ1: S => P; k*S*x\nS = 1; P = 0; k = 1");
    assert_not!(r, RU_RULE_LOOP);
}

#[test]
fn rule_and_event_set_the_same_variable() {
    let r = run("x := 2*S\nJ1: S => P; k*S*x\nE1: at time > 5: x = 3\nS = 1; P = 0; k = 1");
    assert_has!(r, RU_SET_TWICE);
}

#[test]
fn species_in_reactions_set_by_a_rule() {
    let mut d = doc("J1: S => P; k*S\nS = 1; P = 0; k = 1");
    d.model.as_mut().unwrap().rules.push(sbml_io::Rule::Rate {
        variable: "P".into(),
        math: Some(sbml_io::Math::real(1.0)),
        meta: Default::default(),
    });
    let r = mc_core::check(&sbml_io::write(&d));
    assert_has!(r, RU_SET_TWICE);
}

#[test]
fn constant_changed_by_an_event() {
    let mut d = doc("J1: S => P; k*S\nE1: at time > 5: k = 2\nS = 1; P = 0; k = 1");
    let m = d.model.as_mut().unwrap();
    m.parameters.iter_mut().find(|p| p.id == "k").unwrap().constant = Some(true);
    let r = mc_core::check(&sbml_io::write(&d));
    assert_has!(r, RU_CONSTANT_CHANGED);
}

#[test]
fn trigger_true_at_start_and_never_firing() {
    let r = run("J1: S => P; k*S\nE1: at S > 0.5: k = 2\nS = 1; P = 0; k = 1");
    assert_has!(r, RU_TRIGGER_TRUE_AT_START);
    let r = run("J1: S => P; k*S\nE1: at S < 0.5: k = 2\nS = 1; P = 0; k = 1");
    assert_not!(r, RU_TRIGGER_TRUE_AT_START);
    assert_not!(r, RU_EVENT_NEVER);
    // q never changes and q > 5 is false: the event can never fire.
    let r = run("J1: S => P; k*S\nE1: at q > 5: k = 2\nS = 1; P = 0; k = 1; q = 1");
    assert_has!(r, RU_EVENT_NEVER);
}

// ------------------------------------------------- annotations

/// An RDF annotation with one `bqbiol:<qualifier>` link.
fn link(metaid: &str, qualifier: &str, uri: &str) -> String {
    format!(
        "<rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\" \
         xmlns:bqbiol=\"http://biomodels.net/biology-qualifiers/\">\
         <rdf:Description rdf:about=\"#{metaid}\"><bqbiol:{qualifier}><rdf:Bag>\
         <rdf:li rdf:resource=\"{uri}\"/></rdf:Bag></bqbiol:{qualifier}></rdf:Description></rdf:RDF>"
    )
}

/// The simple S => P model with S annotated as `uri` (qualifier `q`).
fn run_annotated(q: &str, uri: &str) -> Report {
    let mut d = doc("J1: S => P; k*S\nS = 1; P = 0; k = 1");
    let s = &mut d.model.as_mut().unwrap().species[0];
    s.meta.metaid = Some("meta_S".into());
    s.meta.annotation = Some(link("meta_S", q, uri));
    mc_core::check(&sbml_io::write(&d))
}

fn message_of<'r>(r: &'r Report, code: &str) -> &'r str {
    r.findings
        .iter()
        .find(|f| f.code == code)
        .map(|f| f.message.as_str())
        .unwrap_or_else(|| panic!("no {code} in:\n{}", show(r)))
}

#[test]
fn missing_database_links() {
    let r = run_annotated("is", "https://identifiers.org/CHEBI:15422");
    let m = message_of(&r, MD_NO_LINKS);
    assert!(m.starts_with("1 of 2 species have no database links"), "{m}");
    assert!(m.contains(": P."), "{m}");
    assert_has!(r, MD_MODEL_INFO);
}

#[test]
fn good_links_pass() {
    for uri in [
        "https://identifiers.org/CHEBI:15422",
        "http://identifiers.org/chebi/CHEBI:15422",
        "https://identifiers.org/kegg.compound:C00002",
    ] {
        let r = run_annotated("is", uri);
        assert_not!(r, MD_BROKEN_LINK);
        assert_not!(r, MD_WRONG_KIND);
    }
}

#[test]
fn broken_links() {
    // ChEBI ids carry their prefix: 15422 alone is not one.
    let r = run_annotated("is", "http://identifiers.org/chebi/15422");
    assert!(message_of(&r, MD_BROKEN_LINK).contains("is not a ChEBI id"), "{}", show(&r));
    let r = run_annotated("is", "https://identifiers.org/refseq_synonym/b0001");
    assert!(message_of(&r, MD_BROKEN_LINK).contains("refseq_synonym (1)"), "{}", show(&r));
}

#[test]
fn old_and_other_links() {
    let r = run_annotated("is", "urn:miriam:chebi:CHEBI%3A15422");
    assert!(message_of(&r, MD_OLD_LINK).contains("urn:miriam"), "{}", show(&r));
    assert_not!(r, MD_BROKEN_LINK);
    let r = run_annotated("is", "http://example.org/glucose");
    assert_has!(r, MD_OTHER_LINK);
}

#[test]
fn link_of_the_wrong_kind() {
    // A species cannot be a reaction.
    let r = run_annotated("is", "https://identifiers.org/rhea:10000");
    assert_has!(r, MD_WRONG_KIND);
    // But it may be part of one, or a version of something.
    let r = run_annotated("isPartOf", "https://identifiers.org/rhea:10000");
    assert_not!(r, MD_WRONG_KIND);
}

#[test]
fn elements_without_names() {
    let r = run("J1: species_0000001 => species_0000002; k*species_0000001\n\
                 species_0000001 = 1; species_0000002 = 0; k = 1");
    assert_has!(r, MD_NO_NAME);
    let r = run("J1: glucose => g6p; k*glucose\nglucose = 1; g6p = 0; k = 1");
    assert_not!(r, MD_NO_NAME);
}
