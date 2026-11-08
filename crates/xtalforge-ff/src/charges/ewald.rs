//! Periodic Coulomb energy and forces by Ewald summation (Task 4.3.5).
//!
//! Ported from Crystal-Builder `xtal/ff/ewald.py`.  Four pieces are all
//! needed for the answer to be a property of the crystal rather than of the
//! cutoff:
//!
//! * real space — `q_i q_j erfc(alpha r) / r` over images inside the cutoff;
//! * reciprocal space — the structure factor on every reciprocal vector
//!   inside `k_cutoff`;
//! * self energy — each charge against its own screening cloud;
//! * the neutralising background — a net-charged cell would otherwise be
//!   infinitely repulsive, so a uniform background is assumed and its energy
//!   removed.
//!
//! `alpha` follows from the real-space cutoff and the requested accuracy, and
//! the reciprocal cutoff follows from `alpha`; halving the real cutoff gives a
//! wider screening and more k-vectors, and the same energy.

use serde::{Deserialize, Serialize};
use xtalforge_core::{Lattice, Structure};

use crate::Error;
use crate::charges::{COULOMB_EV_ANGSTROM, periodic::image_range};
use crate::math::{add, dot, norm, sub, vec_mat};

/// Options for one Ewald evaluation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EwaldOptions {
    /// Real-space cutoff in Angstrom.
    pub cutoff: f64,
    /// Relative accuracy of the split, in `(0, 1)`.
    pub accuracy: f64,
    /// Dielectric screening of the Coulomb interaction.
    pub dielectric: f64,
}

impl Default for EwaldOptions {
    fn default() -> Self {
        Self {
            cutoff: 12.0,
            accuracy: 1e-6,
            dielectric: 1.0,
        }
    }
}

impl EwaldOptions {
    /// Checks that the options are usable.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a non-positive cutoff, an accuracy
    /// outside `(0, 1)`, or a non-positive dielectric.
    pub fn validate(&self) -> Result<(), Error> {
        if !self.cutoff.is_finite() || self.cutoff <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "Ewald cutoff must be positive and finite, got {}",
                self.cutoff
            )));
        }
        if !(self.accuracy.is_finite() && 0.0 < self.accuracy && self.accuracy < 1.0) {
            return Err(Error::InvalidInput(format!(
                "Ewald accuracy must be in (0, 1), got {}",
                self.accuracy
            )));
        }
        if !self.dielectric.is_finite() || self.dielectric <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "dielectric must be positive and finite, got {}",
                self.dielectric
            )));
        }
        Ok(())
    }
}

/// Energy and forces from one Ewald evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct EwaldResult {
    /// Coulomb energy in eV.
    pub energy: f64,
    /// Force on each P1 atom in eV/Angstrom.
    pub forces: Vec<[f64; 3]>,
}

/// Coulomb energy and forces for point charges in a periodic structure.
///
/// `charges` has one entry per P1 atom, in [`xtalforge_core::P1Cell::atoms`]
/// order; for a P1 structure that is the site order.  Dummy atoms participate
/// like any other point charge, so pass zero for them.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for bad options, a charge count that does
/// not match the P1 cell, or a cell with no volume, and a core error when the
/// P1 expansion fails.
pub fn ewald_energy_forces(
    structure: &Structure,
    charges: &[f64],
    options: &EwaldOptions,
) -> Result<EwaldResult, Error> {
    options.validate()?;
    let cell = structure.p1()?;
    if charges.len() != cell.n_atoms() {
        return Err(Error::InvalidInput(format!(
            "got {} charges for {} P1 atoms",
            charges.len(),
            cell.n_atoms()
        )));
    }
    let setup = Setup::new(&cell.lattice, options.cutoff, options.accuracy)?;
    let positions = cell.lattice.to_cart_many(&cell.frac());
    let prefactor = COULOMB_EV_ANGSTROM / options.dielectric;
    let mut forces = vec![[0.0; 3]; positions.len()];
    let mut energy = real_space(
        &cell.lattice,
        &positions,
        charges,
        &setup,
        prefactor,
        &mut forces,
    );
    energy += reciprocal(&positions, charges, &setup, prefactor, &mut forces);
    energy += self_and_background(charges, &setup, prefactor);
    Ok(EwaldResult { energy, forces })
}

/// Convenience wrapper returning only the energy.
///
/// # Errors
///
/// Returns the errors of [`ewald_energy_forces`].
pub fn ewald_energy(
    structure: &Structure,
    charges: &[f64],
    options: &EwaldOptions,
) -> Result<f64, Error> {
    Ok(ewald_energy_forces(structure, charges, options)?.energy)
}

/// The split chosen once for a cell.
#[derive(Debug, Clone)]
pub(crate) struct Setup {
    /// Screening width in inverse Angstrom.
    pub alpha: f64,
    /// Real-space cutoff in Angstrom.
    pub cutoff: f64,
    /// Reciprocal vectors in inverse Angstrom, with `2 pi`.
    pub k_vectors: Vec<[f64; 3]>,
    /// `exp(-k^2 / 4 alpha^2) / k^2` for each vector.
    pub k_factors: Vec<f64>,
    /// Cell volume in Angstrom cubed.
    pub volume: f64,
}

impl Setup {
    pub(crate) fn new(lattice: &Lattice, cutoff: f64, accuracy: f64) -> Result<Self, Error> {
        let volume = lattice.volume().abs();
        if volume <= 0.0 {
            return Err(Error::InvalidInput(
                "Ewald summation needs a cell with positive volume".to_owned(),
            ));
        }
        let tolerance = (-accuracy.ln()).sqrt();
        let alpha = tolerance / cutoff;
        let k_cutoff = 2.0 * alpha * tolerance;
        let reciprocal = lattice.reciprocal(true)?;
        let matrix = reciprocal.matrix();
        let lengths = [norm(matrix[0]), norm(matrix[1]), norm(matrix[2])];
        let mut limits = [1_i32; 3];
        for (index, length) in lengths.iter().enumerate() {
            let value = (k_cutoff / length.max(1e-12)).ceil().max(1.0);
            if let Some(slot) = limits.get_mut(index) {
                *slot = value as i32;
            }
        }
        let mut k_vectors = Vec::new();
        let mut k_factors = Vec::new();
        for h in -limits[0]..=limits[0] {
            for k in -limits[1]..=limits[1] {
                for l in -limits[2]..=limits[2] {
                    let vector = vec_mat([f64::from(h), f64::from(k), f64::from(l)], matrix);
                    let k2 = dot(vector, vector);
                    if k2 > 1e-12 && k2 < k_cutoff * k_cutoff {
                        k_vectors.push(vector);
                        k_factors.push((-k2 / (4.0 * alpha * alpha)).exp() / k2);
                    }
                }
            }
        }
        Ok(Self {
            alpha,
            cutoff,
            k_vectors,
            k_factors,
            volume,
        })
    }
}

/// Real-space half of the Ewald sum, adding forces as it goes.
fn real_space(
    lattice: &Lattice,
    positions: &[[f64; 3]],
    charges: &[f64],
    setup: &Setup,
    prefactor: f64,
    forces: &mut [[f64; 3]],
) -> f64 {
    let n = positions.len();
    let range = image_range(lattice, setup.cutoff);
    let matrix = lattice.matrix();
    let mut energy = 0.0;
    for na in -range[0]..=range[0] {
        for nb in -range[1]..=range[1] {
            for nc in -range[2]..=range[2] {
                let offset = vec_mat([f64::from(na), f64::from(nb), f64::from(nc)], matrix);
                let same_cell = na == 0 && nb == 0 && nc == 0;
                for i in 0..n {
                    for j in 0..n {
                        if same_cell && i == j {
                            continue;
                        }
                        if j < i {
                            continue;
                        }
                        let (Some(first), Some(second)) = (positions.get(i), positions.get(j))
                        else {
                            continue;
                        };
                        let delta = sub(add(*second, offset), *first);
                        let distance = norm(delta);
                        if distance <= 1e-9 || distance > setup.cutoff {
                            continue;
                        }
                        let (Some(qi), Some(qj)) = (charges.get(i), charges.get(j)) else {
                            continue;
                        };
                        let pair = prefactor * qi * qj;
                        let screened = libm_erfc(setup.alpha * distance);
                        // The double lattice sum counts an atom's images on
                        // both sides; an i == j term therefore carries half.
                        let weight = if i == j { 0.5 } else { 1.0 };
                        energy += weight * pair * screened / distance;

                        let gaussian = 2.0 * setup.alpha / std::f64::consts::PI.sqrt()
                            * (-(setup.alpha * distance).powi(2)).exp();
                        let derivative = -(screened / (distance * distance) + gaussian / distance);
                        // `derivative` is dE/dr; the force is its negative.
                        let scale = -pair * derivative / distance;
                        let force = crate::math::scale(delta, scale);
                        if let Some(slot) = forces.get_mut(j) {
                            *slot = add(*slot, force);
                        }
                        if let Some(slot) = forces.get_mut(i) {
                            *slot = sub(*slot, force);
                        }
                    }
                }
            }
        }
    }
    energy
}

/// Reciprocal-space half of the Ewald sum, adding forces as it goes.
fn reciprocal(
    positions: &[[f64; 3]],
    charges: &[f64],
    setup: &Setup,
    prefactor: f64,
    forces: &mut [[f64; 3]],
) -> f64 {
    let n = positions.len();
    if setup.k_vectors.is_empty() || n == 0 {
        return 0.0;
    }
    let mut energy = 0.0;
    for (vector, factor) in setup.k_vectors.iter().zip(&setup.k_factors) {
        let weight = 2.0 * std::f64::consts::PI / setup.volume * factor;
        let mut real = 0.0;
        let mut imaginary = 0.0;
        for (position, charge) in positions.iter().zip(charges) {
            let phase = dot(*vector, *position);
            real += charge * phase.cos();
            imaginary += charge * phase.sin();
        }
        let magnitude = real * real + imaginary * imaginary;
        energy += prefactor * weight * magnitude;

        for (index, (position, charge)) in positions.iter().zip(charges).enumerate() {
            let phase = dot(*vector, *position);
            let sine = phase.sin();
            let cosine = phase.cos();
            // Im(exp(i k.r) conj(S)) = sin(phase) real - cos(phase) imaginary.
            let inner = sine * real - cosine * imaginary;
            let coefficient = -2.0 * charge * inner * weight;
            if let Some(slot) = forces.get_mut(index) {
                let contribution = crate::math::scale(*vector, -prefactor * coefficient);
                *slot = add(*slot, contribution);
            }
        }
    }
    energy
}

/// Self-energy and neutralising background, neither of which moves an atom.
fn self_and_background(charges: &[f64], setup: &Setup, prefactor: f64) -> f64 {
    let sum_squares: f64 = charges.iter().map(|charge| charge * charge).sum();
    let net: f64 = charges.iter().sum();
    let self_energy = -setup.alpha / std::f64::consts::PI.sqrt() * sum_squares;
    let background =
        -std::f64::consts::PI * net * net / (2.0 * setup.alpha * setup.alpha * setup.volume);
    prefactor * (self_energy + background)
}

/// Periodic `1/r` between every pair of atoms, screening cloud included.
///
/// The matrix Ewald and EQeq share: with charges `q`, `E = prefactor / 2 *
/// q^T M q` is [`ewald_energy_forces`] without exclusions.  The diagonal is an
/// atom's interaction with its own images, and the background term is present
/// in every entry as a constant — for a neutral charge set it adds nothing,
/// and what it buys is that every entry is independent of `alpha`.
pub(crate) fn pair_matrix(
    lattice: &Lattice,
    positions: &[[f64; 3]],
    cutoff: f64,
    accuracy: f64,
) -> Result<Vec<Vec<f64>>, Error> {
    let setup = Setup::new(lattice, cutoff, accuracy)?;
    let n = positions.len();
    let mut result = vec![vec![0.0; n]; n];
    let range = image_range(lattice, cutoff);
    let matrix = lattice.matrix();
    for na in -range[0]..=range[0] {
        for nb in -range[1]..=range[1] {
            for nc in -range[2]..=range[2] {
                let offset = vec_mat([f64::from(na), f64::from(nb), f64::from(nc)], matrix);
                for i in 0..n {
                    for j in 0..n {
                        let (Some(first), Some(second)) = (positions.get(i), positions.get(j))
                        else {
                            continue;
                        };
                        let distance = norm(sub(add(*second, offset), *first));
                        if distance <= 1e-9 || distance > cutoff {
                            continue;
                        }
                        let value = libm_erfc(setup.alpha * distance) / distance;
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

    if !setup.k_vectors.is_empty() {
        let weight: Vec<f64> = setup
            .k_factors
            .iter()
            .map(|factor| 4.0 * std::f64::consts::PI / setup.volume * factor)
            .collect();
        for i in 0..n {
            for j in 0..n {
                let (Some(first), Some(second)) = (positions.get(i), positions.get(j)) else {
                    continue;
                };
                let mut value = 0.0;
                for (vector, weight) in setup.k_vectors.iter().zip(&weight) {
                    value += dot(*vector, sub(*first, *second)).cos() * weight;
                }
                if let Some(row) = result.get_mut(i)
                    && let Some(slot) = row.get_mut(j)
                {
                    *slot += value;
                }
            }
        }
    }

    let self_term = 2.0 * setup.alpha / std::f64::consts::PI.sqrt();
    let background = std::f64::consts::PI / (setup.alpha * setup.alpha * setup.volume);
    for (i, row) in result.iter_mut().enumerate() {
        for (j, slot) in row.iter_mut().enumerate() {
            *slot -= background;
            if i == j {
                *slot -= self_term;
            }
        }
    }
    Ok(result)
}

/// Complementary error function without an external special-function crate.
///
/// `erfc` is evaluated through the complementary error function identity and
/// a Chebyshev-free rational approximation is avoided on purpose: the values
/// that matter here are `alpha r` between zero and about six, where the
/// classic Numerical Recipes `erfc` rational form is accurate to `1.2e-7`
/// absolute.  That is two orders below the default Ewald accuracy after the
/// prefactor, and it is deterministic.
pub(crate) fn libm_erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let polynomial = -z * z - 1.265_512_23
        + t * (1.000_023_68
            + t * (0.374_091_96
                + t * (0.096_784_18
                    + t * (-0.186_288_06
                        + t * (0.278_868_07
                            + t * (-1.135_203_98
                                + t * (1.488_515_87 + t * (-0.822_152_23 + t * 0.170_872_77))))))));
    let ans = t * polynomial.exp();
    if x >= 0.0 { ans } else { 2.0 - ans }
}
