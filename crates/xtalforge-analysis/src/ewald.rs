#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Ewald electrostatics for defect and ordering energies (Task 9.5.2).
//!
//! The implementation is the standard real-space plus reciprocal-space sum
//! with a Gaussian screening factor, in eV Angstrom units. It is used to rank
//! orderings and is not a replacement for the full force-field Ewald module.

use xtalforge_core::{Error, Structure, neighbors};

const COULOMB: f64 = 14.399_645_478_6;

/// Options for [`ewald_energy`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct EwaldOptions {
    /// Relative accuracy target for the splitting parameter.
    pub accuracy: f64,
    /// Minimum real-space cutoff in Angstrom.
    pub real_cutoff: f64,
}

impl Default for EwaldOptions {
    fn default() -> Self {
        Self {
            accuracy: 1e-8,
            real_cutoff: 8.0,
        }
    }
}

/// Total electrostatic energy in eV.
///
/// The cell must be charge neutral; the constant background term of a charged
/// cell is deliberately not included.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-neutral cell, a charge count that
/// does not match the P1 atoms, or invalid options.
pub fn ewald_energy(
    structure: &Structure,
    charges: &[f64],
    options: &EwaldOptions,
) -> Result<f64, Error> {
    if !options.real_cutoff.is_finite() || options.real_cutoff <= 0.0 {
        return Err(Error::InvalidInput(
            "the Ewald real-space cutoff must be positive and finite".to_owned(),
        ));
    }
    if !options.accuracy.is_finite() || !(0.0..1.0).contains(&options.accuracy) {
        return Err(Error::InvalidInput(
            "the Ewald accuracy must be in (0, 1)".to_owned(),
        ));
    }
    let p1 = structure.p1()?;
    if charges.len() != p1.atoms.len() {
        return Err(Error::InvalidInput(format!(
            "Ewald needs one charge per atom, got {} charges for {} atoms",
            charges.len(),
            p1.atoms.len()
        )));
    }
    let total: f64 = charges.iter().sum();
    if total.abs() > 1e-4 {
        return Err(Error::InvalidInput(format!(
            "Ewald energy needs a neutral cell, total charge is {total}"
        )));
    }
    if !charges.iter().all(|charge| charge.is_finite()) {
        return Err(Error::InvalidInput(
            "Ewald charges must be finite".to_owned(),
        ));
    }
    let volume = p1.lattice.volume();
    let lengths = p1.lattice.lengths();
    let minimum_width = lengths
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min)
        .max(0.5);
    let alpha = 3.2 / options.real_cutoff.max(2.0 * minimum_width);
    let reciprocal = p1.lattice.reciprocal(false)?;
    let reciprocal_matrix = reciprocal.matrix();
    let k_cut = 2.0 * alpha * (1.0 / options.accuracy).ln().sqrt();
    let max_length = lengths.iter().copied().fold(0.0_f64, f64::max);
    let limit =
        ((k_cut * max_length / (2.0 * std::f64::consts::PI)).ceil() as i32 + 1).clamp(1, 24);
    let mut energy = 0.0;
    // Real space.
    let snapshot = structure.p1()?;
    let mut builder = xtalforge_core::StructureBuilder::new("p1", snapshot.lattice);
    for atom in &snapshot.atoms {
        builder = builder.site(
            xtalforge_core::Site::new(atom.label.clone(), &atom.element, atom.frac)?
                .with_occupancy(atom.occupancy)?,
        );
    }
    let real_structure = builder.build()?;
    for pair in neighbors(&real_structure, options.real_cutoff)? {
        let (Some(&first), Some(&second)) = (charges.get(pair.i), charges.get(pair.j)) else {
            continue;
        };
        let value = first * second * erfc(alpha * pair.distance) / pair.distance;
        if pair.i == pair.j {
            energy += 0.5 * value;
        } else {
            energy += value;
        }
    }
    // Reciprocal space.
    let two_pi = 2.0 * std::f64::consts::PI;
    for h in -limit..=limit {
        for k in -limit..=limit {
            for l in -limit..=limit {
                if h == 0 && k == 0 && l == 0 {
                    continue;
                }
                let vector = [
                    two_pi
                        * (f64::from(h) * reciprocal_matrix[0][0]
                            + f64::from(k) * reciprocal_matrix[1][0]
                            + f64::from(l) * reciprocal_matrix[2][0]),
                    two_pi
                        * (f64::from(h) * reciprocal_matrix[0][1]
                            + f64::from(k) * reciprocal_matrix[1][1]
                            + f64::from(l) * reciprocal_matrix[2][1]),
                    two_pi
                        * (f64::from(h) * reciprocal_matrix[0][2]
                            + f64::from(k) * reciprocal_matrix[1][2]
                            + f64::from(l) * reciprocal_matrix[2][2]),
                ];
                let norm =
                    (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]).sqrt();
                if norm > k_cut || norm < 1e-12 {
                    continue;
                }
                let mut real = 0.0;
                let mut imaginary = 0.0;
                for (atom, &charge) in snapshot.atoms.iter().zip(charges) {
                    let cart = snapshot.lattice.to_cart(atom.frac);
                    let phase = vector[0] * cart[0] + vector[1] * cart[1] + vector[2] * cart[2];
                    real += charge * phase.cos();
                    imaginary += charge * phase.sin();
                }
                energy += 2.0 * std::f64::consts::PI / volume
                    * (-(norm * norm) / (4.0 * alpha * alpha)).exp()
                    / (norm * norm)
                    * (real * real + imaginary * imaginary);
            }
        }
    }
    energy *= COULOMB;
    // Self-energy correction.
    let self_energy: f64 = charges.iter().map(|charge| charge * charge).sum::<f64>() * alpha
        / std::f64::consts::PI.sqrt();
    Ok(energy - COULOMB * self_energy)
}

/// Error-function complement via a Chebyshev-like rational approximation.
fn erfc(value: f64) -> f64 {
    if value < 0.0 {
        return 2.0 - erfc(-value);
    }
    if value > 27.0 {
        return 0.0;
    }
    let t = 1.0 / (1.0 + 0.327_591_1 * value);
    let poly = -value * value - 1.265_512_23
        + t * (1.000_023_68
            + t * (0.374_091_96
                + t * (0.096_784_18
                    + t * (-0.186_288_06
                        + t * (0.278_868_07
                            + t * (-1.135_203_98
                                + t * (1.488_515_87 + t * (-0.822_152_23 + t * 0.170_872_77))))))));
    t * poly.exp()
}
