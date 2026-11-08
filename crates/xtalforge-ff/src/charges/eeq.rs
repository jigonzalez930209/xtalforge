//! EEQ-style charge equilibration (Task 4.3.3).
//!
//! EEQ is electronegativity equilibration with an explicit Coulomb radius per
//! element: the diagonal of the interaction matrix is the chemical hardness and
//! the off-diagonal is the Ohno-Klopman shielded Coulomb kernel
//!
//! ```text
//! E(Q) = sum_i (chi_i Q_i + 1/2 J_i Q_i^2)
//!      + sum_{i<j} C / (eps sqrt(r_ij^2 + R_i R_j)) Q_i Q_j,
//! ```
//!
//! minimized with one Lagrange multiplier for the requested total charge.
//! `eps` is [`ChargeOptions::dielectric`], which QEq ignores and EEQ honours.
//!
//! **Documented approximation.**  The fitted EEQ2019 `chi`/`eta`/`rad` table
//! shipped by sci-form (`src/charges_eeq/eeq2019_data.rs`, Apache-2.0) is *not*
//! vendored in this crate.  This is the documented EEQ kernel with the local
//! table instead: `chi` and `J` come from the existing QEq table in
//! [`crate::charges::data`], and each Coulomb radius is the Ohno radius implied
//! by that hardness, `R_i = e^2 / (4 pi eps0 J_i)`.  Deriving the radius from
//! the same hardness keeps the kernel's `r = 0` limit below the diagonal and
//! the linear system well conditioned; an independently fitted radius without
//! refitting the hardness would over-screen and can make water singular.
//! Because QEq's shield is the same Ohno radius, the default options
//! (`dielectric = 1`) reproduce [`crate::charges::qeq_charges`] to rounding;
//! raising the dielectric through [`ChargeOptions::eqeq`] is the documented way
//! to obtain the screened EEQ-style charges.
//!
//! The real-space sum runs over periodic images inside `options.cutoff`, exactly
//! like QEq, so a small cell sees every copy of a neighbour and the result is a
//! property of the crystal rather than of the cell's orientation.

use xtalforge_core::{Lattice, Structure, element_info};

use crate::Error;
use crate::charges::data::qeq_params;
use crate::charges::periodic::image_range;
use crate::charges::qeq::solve_constrained;
use crate::charges::{COULOMB_EV_ANGSTROM, ChargeOptions, Charges, MAX_CHARGE_ATOMS};
use crate::math::{add, norm, sub, vec_mat};

/// One element's EEQ parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
struct EeqParams {
    /// Electronegativity in eV.
    chi: f64,
    /// Chemical hardness in eV.
    hardness: f64,
    /// Ohno Coulomb radius in Angstrom.
    radius: f64,
}

/// Reads the local table and derives the Coulomb radius.
fn eeq_params(element: &str) -> Result<EeqParams, Error> {
    let params = qeq_params(element)?;
    if !params.hardness.is_finite() || params.hardness <= 0.0 {
        return Err(Error::Data(format!(
            "EEQ hardness for {element} is not positive and finite"
        )));
    }
    Ok(EeqParams {
        chi: params.chi,
        hardness: params.hardness,
        radius: COULOMB_EV_ANGSTROM / params.hardness,
    })
}

/// Solves EEQ for one structure.
pub(crate) fn solve(structure: &Structure, options: &ChargeOptions) -> Result<Charges, Error> {
    let cell = structure.p1()?;
    let n = cell.n_atoms();
    if n == 0 {
        return Ok(Charges::empty());
    }
    if n > MAX_CHARGE_ATOMS {
        return Err(Error::InvalidInput(format!(
            "EEQ solves a dense {n} x {n} system; above {MAX_CHARGE_ATOMS} atoms set site charges instead"
        )));
    }

    let mut indices = Vec::new();
    let mut chi = Vec::new();
    let mut hardness = Vec::new();
    let mut radius = Vec::new();
    let mut positions = Vec::new();
    for (index, atom) in cell.atoms.iter().enumerate() {
        if element_info(&atom.element)?.is_dummy() {
            continue;
        }
        let params = eeq_params(&atom.element)?;
        indices.push(index);
        chi.push(params.chi);
        hardness.push(params.hardness);
        radius.push(params.radius);
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

    let interaction = interaction_matrix(
        &cell.lattice,
        &positions,
        &hardness,
        &radius,
        options.cutoff,
        options.dielectric,
    );
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

/// The EEQ interaction matrix in eV per unit charge product.
///
/// The diagonal is the hardness; every off-diagonal is the Ohno-Klopman kernel
/// `C / (eps sqrt(r^2 + R_i R_j))` summed over the periodic images inside the
/// cutoff.  The `eps` scaling is what separates EEQ's screening from QEq's
/// unscreened Ohno kernel.
fn interaction_matrix(
    lattice: &Lattice,
    positions: &[[f64; 3]],
    hardness: &[f64],
    radius: &[f64],
    cutoff: f64,
    dielectric: f64,
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
    let prefactor = COULOMB_EV_ANGSTROM / dielectric;
    for na in -range[0]..=range[0] {
        for nb in -range[1]..=range[1] {
            for nc in -range[2]..=range[2] {
                let offset = vec_mat([f64::from(na), f64::from(nb), f64::from(nc)], base);
                let same_cell = na == 0 && nb == 0 && nc == 0;
                for i in 0..n {
                    for j in 0..n {
                        if same_cell && i == j {
                            continue;
                        }
                        let (Some(first), Some(second)) = (positions.get(i), positions.get(j))
                        else {
                            continue;
                        };
                        let distance = norm(sub(add(*second, offset), *first));
                        if distance > cutoff {
                            continue;
                        }
                        let (Some(ri), Some(rj)) = (radius.get(i), radius.get(j)) else {
                            continue;
                        };
                        let width = (ri * rj).sqrt();
                        let value = prefactor / (distance * distance + width * width).sqrt();
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
