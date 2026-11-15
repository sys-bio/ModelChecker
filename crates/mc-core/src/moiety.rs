//! Conserved moieties, used to rule out states a simulation cannot reach.
//!
//! If `c >= 0`, `c != 0`, and `c^T N = 0` for the stoichiometry matrix `N`,
//! then `c . x` is constant. When its initial value is positive, the species
//! in `c`'s support can never all be zero at once. A candidate witness that
//! sets them all to zero is therefore not a reachable state.
//!
//! Only species changed by reactions alone take part; a species also set by
//! a rule or event can break any conservation it is in.

use std::cell::RefCell;
use std::collections::HashMap;

pub struct Conservation {
    /// Species that only reactions change, with their initial values.
    index: HashMap<String, usize>,
    x0: Vec<f64>,
    /// Columns of N: per reaction, (species index, net stoichiometry).
    columns: Vec<Vec<(usize, f64)>>,
    cache: RefCell<HashMap<Vec<usize>, bool>>,
}

const TOL: f64 = 1e-9;

impl Conservation {
    pub fn new(species: Vec<(String, f64)>, reactions: Vec<Vec<(String, f64)>>) -> Self {
        let index: HashMap<String, usize> = species
            .iter()
            .enumerate()
            .map(|(i, (s, _))| (s.clone(), i))
            .collect();
        let x0 = species.iter().map(|(_, v)| *v).collect();
        let columns = reactions
            .into_iter()
            .map(|col| {
                let mut merged: Vec<(usize, f64)> = Vec::new();
                for (s, n) in col {
                    let Some(&i) = index.get(&s) else { continue };
                    match merged.iter_mut().find(|(j, _)| *j == i) {
                        Some(e) => e.1 += n,
                        None => merged.push((i, n)),
                    }
                }
                merged.retain(|(_, n)| n.abs() > TOL);
                merged
            })
            .filter(|c| !c.is_empty())
            .collect();
        Conservation {
            index,
            x0,
            columns,
            cache: RefCell::new(HashMap::new()),
        }
    }

    /// Whether the species in `zeros` cannot all be zero at once, because a
    /// conserved total with a positive value lies within them.
    pub fn cannot_all_be_zero<'a>(&self, zeros: impl IntoIterator<Item = &'a str>) -> bool {
        let mut z: Vec<usize> = zeros
            .into_iter()
            .filter_map(|s| self.index.get(s).copied())
            .collect();
        z.sort_unstable();
        z.dedup();
        if z.is_empty() {
            return false;
        }
        if let Some(&r) = self.cache.borrow().get(&z) {
            return r;
        }
        let r = self.compute(&z);
        self.cache.borrow_mut().insert(z, r);
        r
    }

    fn compute(&self, z: &[usize]) -> bool {
        // A: one row per reaction, one column per species in z, so that
        // A c = 0 is c^T N_z = 0.
        let k = z.len();
        let pos: HashMap<usize, usize> = z.iter().enumerate().map(|(c, &s)| (s, c)).collect();
        let mut a: Vec<Vec<f64>> = self
            .columns
            .iter()
            .filter_map(|col| {
                let mut row = vec![0.0; k];
                let mut any = false;
                for &(s, n) in col {
                    if let Some(&c) = pos.get(&s) {
                        row[c] = n;
                        any = true;
                    }
                }
                any.then_some(row)
            })
            .collect();
        // Reduced row echelon form.
        let mut pivots: Vec<(usize, usize)> = Vec::new(); // (row, col)
        let mut row = 0;
        for col in 0..k {
            if row >= a.len() {
                break;
            }
            let (best, val) = (row..a.len())
                .map(|r| (r, a[r][col].abs()))
                .fold((row, 0.0), |b, x| if x.1 > b.1 { x } else { b });
            if val < TOL {
                continue;
            }
            a.swap(row, best);
            let p = a[row][col];
            for v in a[row].iter_mut() {
                *v /= p;
            }
            for r in 0..a.len() {
                if r != row && a[r][col].abs() > TOL {
                    let f = a[r][col];
                    let pivot_row = a[row].clone();
                    for (v, pv) in a[r].iter_mut().zip(&pivot_row) {
                        *v -= f * pv;
                    }
                }
            }
            pivots.push((row, col));
            row += 1;
        }
        let pivot_cols: Vec<usize> = pivots.iter().map(|&(_, c)| c).collect();
        // One null-space basis vector per free column.
        for free in (0..k).filter(|c| !pivot_cols.contains(c)) {
            let mut c = vec![0.0; k];
            c[free] = 1.0;
            for &(r, pc) in &pivots {
                c[pc] = -a[r][free];
            }
            let nonneg = c.iter().all(|v| *v >= -TOL);
            let nonpos = c.iter().all(|v| *v <= TOL);
            if !(nonneg || nonpos) {
                continue;
            }
            let sign = if nonneg { 1.0 } else { -1.0 };
            let total: f64 = c
                .iter()
                .zip(z)
                .map(|(ci, &s)| sign * ci * self.x0[s].max(0.0))
                .sum();
            if total > TOL {
                return true;
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cons(species: &[(&str, f64)], reactions: &[&[(&str, f64)]]) -> Conservation {
        Conservation::new(
            species.iter().map(|(s, v)| (s.to_string(), *v)).collect(),
            reactions
                .iter()
                .map(|r| r.iter().map(|(s, n)| (s.to_string(), *n)).collect())
                .collect(),
        )
    }

    #[test]
    fn cycle_total_cannot_vanish() {
        // E + S <-> ES -> E + P: E + ES is conserved.
        let c = cons(
            &[("E", 1.0), ("S", 10.0), ("ES", 0.0), ("P", 0.0)],
            &[
                &[("E", -1.0), ("S", -1.0), ("ES", 1.0)],
                &[("ES", -1.0), ("E", 1.0), ("P", 1.0)],
            ],
        );
        assert!(c.cannot_all_be_zero(["E", "ES"]));
        assert!(!c.cannot_all_be_zero(["E"]));
        assert!(!c.cannot_all_be_zero(["S", "P"]));
        // S + ES + P is conserved too.
        assert!(c.cannot_all_be_zero(["S", "ES", "P"]));
    }

    #[test]
    fn open_system_can_empty() {
        let c = cons(&[("S", 1.0)], &[&[("S", 1.0)], &[("S", -1.0)]]);
        assert!(!c.cannot_all_be_zero(["S"]));
    }
}
