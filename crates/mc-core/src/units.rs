//! Units as SI base-unit exponents plus a scale, and how a model's unit
//! declarations resolve to them.

use std::collections::HashMap;
use std::fmt;

use sbml_io::Model;

/// metre, kilogram, second, ampere, kelvin, mole, candela, item.
pub const BASE_SYMBOLS: [&str; 8] = ["m", "kg", "s", "A", "K", "mol", "cd", "item"];
const M: usize = 0;
const KG: usize = 1;
const S: usize = 2;
const A: usize = 3;
const K: usize = 4;
const MOL: usize = 5;
const CD: usize = 6;
const ITEM: usize = 7;

/// A unit: `10^log10 * prod(base_i ^ exp_i)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Unit {
    pub exp: [f64; 8],
    pub log10: f64,
}

const TOL: f64 = 1e-9;

impl Unit {
    pub const DIMENSIONLESS: Unit = Unit {
        exp: [0.0; 8],
        log10: 0.0,
    };

    fn base(i: usize, power: f64, log10: f64) -> Unit {
        let mut exp = [0.0; 8];
        exp[i] = power;
        Unit { exp, log10 }
    }

    fn of(parts: &[(usize, f64)], log10: f64) -> Unit {
        let mut exp = [0.0; 8];
        for &(i, p) in parts {
            exp[i] += p;
        }
        Unit { exp, log10 }
    }

    pub fn mul(&self, o: &Unit) -> Unit {
        let mut exp = self.exp;
        for (e, x) in exp.iter_mut().zip(o.exp) {
            *e += x;
        }
        Unit {
            exp,
            log10: self.log10 + o.log10,
        }
    }

    pub fn div(&self, o: &Unit) -> Unit {
        self.mul(&o.pow(-1.0))
    }

    pub fn pow(&self, p: f64) -> Unit {
        let mut exp = self.exp;
        for e in exp.iter_mut() {
            *e *= p;
        }
        Unit {
            exp,
            log10: self.log10 * p,
        }
    }

    /// Same dimensions (exponents), whatever the scale.
    pub fn same_dimensions(&self, o: &Unit) -> bool {
        self.exp.iter().zip(o.exp).all(|(a, b)| (a - b).abs() < TOL)
    }

    /// Same dimensions and the same scale.
    pub fn same(&self, o: &Unit) -> bool {
        self.same_dimensions(o) && (self.log10 - o.log10).abs() < 1e-6
    }

    pub fn is_dimensionless(&self) -> bool {
        self.exp.iter().all(|e| e.abs() < TOL)
    }

    /// The factor between two units of the same dimensions: `self / o`.
    pub fn factor(&self, o: &Unit) -> f64 {
        10f64.powf(self.log10 - o.log10)
    }
}

impl fmt::Display for Unit {
    /// SI text, for example `1000 mol/(m^3*s)`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let term = |i: usize, p: f64| {
            let p = p.abs();
            if (p - 1.0).abs() < TOL {
                BASE_SYMBOLS[i].to_string()
            } else if (p - p.round()).abs() < TOL {
                format!("{}^{}", BASE_SYMBOLS[i], p.round() as i64)
            } else {
                // Fractional powers come from Hill coefficients: four
                // places is plenty, and hides round-off.
                let s = format!("{p:.4}");
                let s = s.trim_end_matches('0').trim_end_matches('.');
                format!("{}^{}", BASE_SYMBOLS[i], s)
            }
        };
        let num: Vec<String> = (0..8)
            .filter(|&i| self.exp[i] > TOL)
            .map(|i| term(i, self.exp[i]))
            .collect();
        let den: Vec<String> = (0..8)
            .filter(|&i| self.exp[i] < -TOL)
            .map(|i| term(i, self.exp[i]))
            .collect();
        let mut s = String::new();
        if self.log10.abs() > 1e-6 {
            let factor = 10f64.powf(self.log10);
            let r = self.log10.round();
            if (self.log10 - r).abs() < 1e-6 {
                s.push_str(&format!("1e{} ", r as i64));
            } else {
                s.push_str(&format!("{} ", crate::sample::fmt_num(factor)));
            }
        }
        s.push_str(&match num.len() {
            0 if den.is_empty() => "dimensionless".to_string(),
            0 => "1".to_string(),
            _ => num.join("*"),
        });
        match den.len() {
            0 => {}
            1 => s.push_str(&format!("/{}", den[0])),
            _ => s.push_str(&format!("/({})", den.join("*"))),
        }
        f.write_str(s.trim())
    }
}

/// The SBML base units.
pub fn base_unit(kind: &str) -> Option<Unit> {
    Some(match kind {
        "dimensionless" | "radian" | "steradian" => Unit::DIMENSIONLESS,
        "avogadro" => Unit {
            exp: [0.0; 8],
            log10: 6.02214076e23f64.log10(),
        },
        "metre" | "meter" => Unit::base(M, 1.0, 0.0),
        "kilogram" => Unit::base(KG, 1.0, 0.0),
        "gram" => Unit::base(KG, 1.0, -3.0),
        "second" => Unit::base(S, 1.0, 0.0),
        "ampere" => Unit::base(A, 1.0, 0.0),
        "kelvin" | "celsius" => Unit::base(K, 1.0, 0.0),
        "mole" => Unit::base(MOL, 1.0, 0.0),
        "candela" => Unit::base(CD, 1.0, 0.0),
        "item" => Unit::base(ITEM, 1.0, 0.0),
        "litre" | "liter" => Unit::base(M, 3.0, -3.0),
        "hertz" | "becquerel" => Unit::base(S, -1.0, 0.0),
        "katal" => Unit::of(&[(MOL, 1.0), (S, -1.0)], 0.0),
        "newton" => Unit::of(&[(KG, 1.0), (M, 1.0), (S, -2.0)], 0.0),
        "pascal" => Unit::of(&[(KG, 1.0), (M, -1.0), (S, -2.0)], 0.0),
        "joule" => Unit::of(&[(KG, 1.0), (M, 2.0), (S, -2.0)], 0.0),
        "watt" => Unit::of(&[(KG, 1.0), (M, 2.0), (S, -3.0)], 0.0),
        "coulomb" => Unit::of(&[(A, 1.0), (S, 1.0)], 0.0),
        "volt" => Unit::of(&[(KG, 1.0), (M, 2.0), (S, -3.0), (A, -1.0)], 0.0),
        "farad" => Unit::of(&[(KG, -1.0), (M, -2.0), (S, 4.0), (A, 2.0)], 0.0),
        "ohm" => Unit::of(&[(KG, 1.0), (M, 2.0), (S, -3.0), (A, -2.0)], 0.0),
        "siemens" => Unit::of(&[(KG, -1.0), (M, -2.0), (S, 3.0), (A, 2.0)], 0.0),
        "weber" => Unit::of(&[(KG, 1.0), (M, 2.0), (S, -2.0), (A, -1.0)], 0.0),
        "tesla" => Unit::of(&[(KG, 1.0), (S, -2.0), (A, -1.0)], 0.0),
        "henry" => Unit::of(&[(KG, 1.0), (M, 2.0), (S, -2.0), (A, -2.0)], 0.0),
        "gray" | "sievert" => Unit::of(&[(M, 2.0), (S, -2.0)], 0.0),
        "lumen" => Unit::base(CD, 1.0, 0.0),
        "lux" => Unit::of(&[(CD, 1.0), (M, -2.0)], 0.0),
        _ => return None,
    })
}

/// A model's unit declarations, resolved.
pub struct Units {
    defs: HashMap<String, Unit>,
    level: u32,
    pub substance: Option<Unit>,
    pub time: Option<Unit>,
    pub volume: Option<Unit>,
    pub area: Option<Unit>,
    pub length: Option<Unit>,
    pub extent: Option<Unit>,
    /// Names for printing: unit definitions by their resolved unit.
    names: Vec<(Unit, String)>,
}

impl Units {
    pub fn new(model: &Model, level: u32) -> Self {
        let mut defs = HashMap::new();
        let mut names = Vec::new();
        for ud in &model.unit_definitions {
            let mut u = Unit::DIMENSIONLESS;
            let mut ok = true;
            for part in &ud.units {
                let Some(b) = base_unit(&part.kind) else {
                    ok = false;
                    break;
                };
                let mult = part.multiplier.unwrap_or(1.0);
                let scale = part.scale.unwrap_or(0) as f64;
                let e = part.exponent.unwrap_or(1.0);
                if mult <= 0.0 {
                    ok = false;
                    break;
                }
                // (multiplier * 10^scale * base)^exponent
                let one = Unit {
                    exp: b.exp,
                    log10: b.log10 + scale + mult.log10(),
                };
                u = u.mul(&one.pow(e));
            }
            if ok {
                defs.insert(ud.id.clone(), u);
                names.push((u, ud.name.clone().unwrap_or_else(|| ud.id.clone())));
            }
        }
        let mut units = Units {
            defs,
            level,
            substance: None,
            time: None,
            volume: None,
            area: None,
            length: None,
            extent: None,
            names,
        };
        let pick = |u: &Units, attr: &Option<String>, l2: &str| -> Option<Unit> {
            if level >= 3 {
                attr.as_deref().and_then(|a| u.resolve(a))
            } else {
                u.resolve(l2)
            }
        };
        units.substance = pick(&units, &model.substance_units, "substance");
        units.time = pick(&units, &model.time_units, "time");
        units.volume = pick(&units, &model.volume_units, "volume");
        units.area = pick(&units, &model.area_units, "area");
        units.length = pick(&units, &model.length_units, "length");
        units.extent = if level >= 3 {
            model.extent_units.as_deref().and_then(|a| units.resolve(a))
        } else {
            units.substance
        };
        units
    }

    /// A unit id: a unit definition, a base unit, or a Level 2 built-in.
    pub fn resolve(&self, id: &str) -> Option<Unit> {
        if let Some(u) = self.defs.get(id) {
            return Some(*u);
        }
        if let Some(u) = base_unit(id) {
            return Some(u);
        }
        if self.level < 3 {
            return match id {
                "substance" => base_unit("mole"),
                "volume" => base_unit("litre"),
                "area" => Some(Unit::base(M, 2.0, 0.0)),
                "length" => base_unit("metre"),
                "time" => base_unit("second"),
                _ => None,
            };
        }
        None
    }

    /// Text for a unit: a unit definition's name if one matches exactly,
    /// else SI.
    pub fn describe(&self, u: &Unit) -> String {
        if let Some((_, n)) = self.names.iter().find(|(d, _)| d.same(u)) {
            return format!("{n} ({u})");
        }
        u.to_string()
    }

    /// Whether the model declares any units at all.
    pub fn any_declared(&self) -> bool {
        !self.defs.is_empty()
            || self.substance.is_some()
            || self.time.is_some()
            || self.extent.is_some()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn litre_and_molar() {
        let l = base_unit("litre").unwrap();
        let mol = base_unit("mole").unwrap();
        let molar = mol.div(&l);
        assert_eq!(molar.to_string(), "1e3 mol/m^3");
        let per_s = base_unit("second").unwrap().pow(-1.0);
        assert_eq!(mol.mul(&per_s).to_string(), "mol/s");
        assert!(molar.mul(&l).same(&mol));
    }
}
