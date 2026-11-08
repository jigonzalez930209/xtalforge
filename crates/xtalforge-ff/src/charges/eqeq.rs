//! EQeq charge equilibration (Task 4.3.2).
//!
//! Ported from Crystal-Builder `xtal/ff/charges/eqeq.py`.  EQeq keeps QEq's
//! energy and changes where `chi` and `J` come from: each element is expanded
//! about the charge it usually carries, and the parabola is read off its
//! measured ionisation energies,
//!
//! ```text
//! chi = (IE_{c+1} + IE_c) / 2 - c J,    J = IE_{c+1} - IE_c,
//! ```
//!
//! with `IE_0` the electron affinity.  The `- c J` moves the expansion back to
//! `Q = 0`; the charges returned are the atoms' charges, not offsets from `c`.
//!
//! **Documented variant.**  The periodic Coulomb matrix is the converged
//! Ewald [`pair_matrix`](crate::charges::ewald) plus the reference's Gaussian
//! overlap correction, and the published dielectric is 1.2
//! ([`ChargeOptions::eqeq`](crate::charges::ChargeOptions::eqeq)).  The table
//! covers the common elements; an element outside it is refused rather than
//! guessed.  The charges are an estimate to look over, not a published result.

use xtalforge_core::{Lattice, Structure, element_info};

use crate::Error;
use crate::charges::data::eqeq_params;
use crate::charges::ewald::pair_matrix;
use crate::charges::periodic::image_range;
use crate::charges::qeq::solve_constrained;
use crate::charges::{COULOMB_EV_ANGSTROM, ChargeOptions, Charges, MAX_CHARGE_ATOMS};
use crate::math::{add, norm, sub, vec_mat};

/// Overlap term is summed until `exp(-(a r)^2)` falls below this.
const OVERLAP_TOLERANCE: f64 = 1e-12;

/// Solves EQeq for one structure.
pub(crate) fn solve(structure: &Structure, options: &ChargeOptions) -> Result<Charges, Error> {
    let cell = structure.p1()?;
    let n = cell.n_atoms();
    if n == 0 {
        return Ok(Charges::empty());
    }
    if n > MAX_CHARGE_ATOMS {
        return Err(Error::InvalidInput(format!(
            "EQeq solves a dense {n} x {n} system; above {MAX_CHARGE_ATOMS} atoms set site charges instead"
        )));
    }

    let mut indices = Vec::new();
    let mut chi = Vec::new();
    let mut hardness = Vec::new();
    let mut positions = Vec::new();
    for (index, atom) in cell.atoms.iter().enumerate() {
        if element_info(&atom.element)?.is_dummy() {
            continue;
        }
        let centre = options.centres.get(&atom.element).copied();
        let (value, j) = eqeq_params(&atom.element, centre)?;
        indices.push(index);
        chi.push(value);
        hardness.push(j);
        positions.push(cell.lattice.to_cart(atom.frac));
    }

    let mut per_atom = vec![0.0; n];
    if indices.is_empty() {
        return Ok(Charges {
            per_atom,
            total: 0.0,
            iterations: 0,
        });
    }

    let coulomb = pair_matrix(&cell.lattice, &positions, options.cutoff, options.accuracy)?;
    let overlap = overlap_matrix(&cell.lattice, &positions, &hardness);
    let scale = options.dielectric * COULOMB_EV_ANGSTROM / 2.0;
    let mut interaction = vec![vec![0.0; positions.len()]; positions.len()];
    for i in 0..positions.len() {
        for j in 0..positions.len() {
            let pair = coulomb
                .get(i)
                .and_then(|row| row.get(j))
                .copied()
                .unwrap_or(0.0);
            let correction = overlap
                .get(i)
                .and_then(|row| row.get(j))
                .copied()
                .unwrap_or(0.0);
            let mut value = scale * (pair + correction);
            if i == j {
                value += hardness.get(i).copied().unwrap_or(0.0);
            }
            if let Some(row) = interaction.get_mut(i)
                && let Some(slot) = row.get_mut(j)
            {
                *slot = value;
            }
        }
    }

    let solution = solve_constrained(&chi, &interaction, options.total_charge)?;
    for (slot, index) in indices.into_iter().enumerate() {
        if let (Some(value), Some(target)) = (solution.get(slot), per_atom.get_mut(index)) {
            *target = *value;
        }
    }
    let total = per_atom.iter().sum();
    Ok(Charges {
        per_atom,
        total,
        iterations: 1,
    })
}

/// The Gaussian overlap correction to `1/r`, summed over periodic images.
///
/// `exp(-(a r)^2) (2a - a^2 r - 1/r)` with `a = sqrt(J_i J_j) / k` cancels the
/// Coulomb `1/r` as two atoms merge and leaves the combined hardness.  An
/// atom's overlap with itself is its hardness, already on the diagonal, so the
/// `r = 0` term is left out.
fn overlap_matrix(lattice: &Lattice, positions: &[[f64; 3]], hardness: &[f64]) -> Vec<Vec<f64>> {
    let n = positions.len();
    let mut result = vec![vec![0.0; n]; n];
    let mut smallest = f64::INFINITY;
    for value in hardness {
        let a = value.sqrt() / COULOMB_EV_ANGSTROM;
        if a > 0.0 {
            smallest = smallest.min(a);
        }
    }
    if !smallest.is_finite() || smallest <= 0.0 {
        return result;
    }
    let reach = (-OVERLAP_TOLERANCE.ln()).sqrt() / smallest;
    let range = image_range(lattice, reach);
    let matrix = lattice.matrix();
    let inverse = lattice.inverse();
    let frac: Vec<[f64; 3]> = positions
        .iter()
        .map(|position| vec_mat(*position, inverse))
        .collect();
    for na in -range[0]..=range[0] {
        for nb in -range[1]..=range[1] {
            for nc in -range[2]..=range[2] {
                let shift = [f64::from(na), f64::from(nb), f64::from(nc)];
                for i in 0..n {
                    for j in 0..n {
                        if na == 0 && nb == 0 && nc == 0 && i == j {
                            continue;
                        }
                        let (Some(first), Some(second)) = (frac.get(i), frac.get(j)) else {
                            continue;
                        };
                        let mut delta = sub(*second, *first);
                        delta = [
                            delta[0] - delta[0].round(),
                            delta[1] - delta[1].round(),
                            delta[2] - delta[2].round(),
                        ];
                        let distance = norm(vec_mat(add(delta, shift), matrix));
                        if distance <= 1e-9 || distance > reach {
                            continue;
                        }
                        let (Some(ji), Some(jj)) = (hardness.get(i), hardness.get(j)) else {
                            continue;
                        };
                        let a = (ji * jj).sqrt() / COULOMB_EV_ANGSTROM;
                        let value = (-(a * distance).powi(2)).exp()
                            * (2.0 * a - a * a * distance - 1.0 / distance);
                        if let Some(row) = result.get_mut(i)
                            && let Some(slot) = row.get_mut(j)
                        {
                            *slot += value;
                        }
                    }
                }
            }
        }
    }
    result
}
