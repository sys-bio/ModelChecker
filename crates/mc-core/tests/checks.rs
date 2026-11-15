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
