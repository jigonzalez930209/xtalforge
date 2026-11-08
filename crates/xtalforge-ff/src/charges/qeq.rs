//! QEq charge equilibration (Task 4.3.1).
//!
//! Ported from Crystal-Builder `xtal/ff/uff/qeq.py`.  The energy is
//!
//! ```text
//! E(Q) = sum_i (chi_i Q_i + 1/2 J_i Q_i^2)
//!      + sum_{i<j} J_ij(r_ij) Q_i Q_j
//! ```
//!
//! and the solution is the set of charges that makes every `dE/dQ_i` equal
//! while the total charge is fixed.  `J_ij` is the Ohno-Klopman shielded
//! form, not the published QEq Slater integral: it tends to the combined
//! hardness as two atoms merge and to `1/r` when they separate, which are the
//! two limits the charges depend on.

use xtalforge_core::{Lattice, Structure, element_info};

use crate::Error;
use crate::charges::data::qeq_params;
use crate::charges::periodic::image_range;
use crate::charges::{COULOMB_EV_ANGSTROM, ChargeOptions, Charges, MAX_CHARGE_ATOMS};
use crate::math::{add, norm, sub, vec_mat};

/// Solves QEq for one structure.
pub(crate) fn solve(structure: &Structure, options: &ChargeOptions) -> Result<Charges, Error> {
    let cell = structure.p1()?;
    let n = cell.n_atoms();
    if n == 0 {
        return Ok(Charges::empty());
    }
    if n > MAX_CHARGE_ATOMS {
        return Err(Error::InvalidInput(format!(
            "QEq solves a dense {n} x {n} system; above {MAX_CHARGE_ATOMS} atoms set site charges instead"
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
        let params = qeq_params(&atom.element)?;
        indices.push(index);
        chi.push(params.chi);
        hardness.push(params.hardness);
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

    let interaction = interaction_matrix(&cell.lattice, &positions, &hardness, options.cutoff);
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

/// The QEq interaction matrix, in eV per unit charge product.
///
/// The diagonal is each atom's hardness.  Off it, the shielded Coulomb term
/// is summed over every image inside `cutoff`, because a small cell sees
/// several copies of the same neighbour and each is a real interaction.
fn interaction_matrix(
    lattice: &Lattice,
    positions: &[[f64; 3]],
    hardness: &[f64],
    cutoff: f64,
) -> Vec<Vec<f64>> {
    let n = positions.len();
    let mut matrix = vec![vec![0.0; n]; n];
    for (index, value) in hardness.iter().enumerate() {
        if let Some(row) = matrix.get_mut(index)
            && let Some(slot) = row.get_mut(index)
        {
            *slot = *value;
        }
    }
    let range = image_range(lattice, cutoff);
    let base = lattice.matrix();
    for na in -range[0]..=range[0] {
        for nb in -range[1]..=range[1] {
            for nc in -range[2]..=range[2] {
                let offset = vec_mat([f64::from(na), f64::from(nb), f64::from(nc)], base);
                let shifted = na == 0 && nb == 0 && nc == 0;
                for i in 0..n {
                    for j in 0..n {
                        if shifted && i == j {
                            continue;
                        }
                        let (Some(first), Some(second)) = (positions.get(i), positions.get(j))
                        else {
                            continue;
                        };
                        let distance = norm(sub(add(*second, offset), *first));
                        if distance <= 1e-9 || distance > cutoff {
                            continue;
                        }
                        let (Some(hi), Some(hj)) = (hardness.get(i), hardness.get(j)) else {
                            continue;
                        };
                        let combined = (hi * hj).sqrt().max(1e-6);
                        let screening = COULOMB_EV_ANGSTROM / combined;
                        let value = COULOMB_EV_ANGSTROM
                            / (distance * distance + screening * screening).sqrt();
                        if let Some(row) = matrix.get_mut(i)
                            && let Some(slot) = row.get_mut(j)
                        {
                            *slot += value;
                        }
                    }
                }
            }
        }
    }
    matrix
}

/// Solves the equality-constrained quadratic program by one KKT system.
pub(crate) fn solve_constrained(
    chi: &[f64],
    interaction: &[Vec<f64>],
    total_charge: f64,
) -> Result<Vec<f64>, Error> {
    let n = chi.len();
    if n == 1 {
        return Ok(vec![total_charge]);
    }
    if interaction.len() != n {
        return Err(Error::Data(
            "QEq interaction matrix has the wrong size".to_owned(),
        ));
    }
    let mut matrix = vec![vec![0.0; n + 1]; n + 1];
    for i in 0..n {
        for j in 0..n {
            let value = interaction
                .get(i)
                .and_then(|row| row.get(j))
                .copied()
                .unwrap_or(0.0);
            if let Some(row) = matrix.get_mut(i) {
                if let Some(slot) = row.get_mut(j) {
                    *slot = value;
                }
                if let Some(slot) = row.get_mut(n) {
                    *slot = 1.0;
                }
            }
        }
        if let Some(row) = matrix.get_mut(n)
            && let Some(slot) = row.get_mut(i)
        {
            *slot = 1.0;
        }
    }
    let mut right = vec![0.0; n + 1];
    for (index, value) in chi.iter().enumerate() {
        if let Some(slot) = right.get_mut(index) {
            *slot = -value;
        }
    }
    if let Some(slot) = right.get_mut(n) {
        *slot = total_charge;
    }
    let solution = crate::math::solve_linear(&mut matrix, &mut right).ok_or_else(|| {
        Error::InvalidInput(
            "QEq linear system is singular; two atoms may share a position".to_owned(),
        )
    })?;
    let mut charges = Vec::with_capacity(n);
    for index in 0..n {
        charges.push(*solution.get(index).unwrap_or(&0.0));
    }
    Ok(charges)
}
