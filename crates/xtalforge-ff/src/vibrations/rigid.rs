//! Mass-weighted translation and rotation vectors and their projector.
//!
//! These are the rigid motions the normal-mode analysis removes from the
//! mass-weighted Hessian.  A crystal contributes three translations; a
//! molecule contributes three translations plus three rotations, of which a
//! linear molecule retains only two because one rotation is the zero vector.
//! The basis is orthonormalised by modified Gram-Schmidt with a fixed order,
//! so the projection is deterministic.

use xtalforge_core::Structure;

use crate::Error;
use crate::math::cross;

use super::hessian::DisplacementMode;

/// Builds the orthonormal mass-weighted rigid-motion basis.
pub(crate) fn rigid_basis(
    structure: &Structure,
    masses: &[f64],
    mode: DisplacementMode,
) -> Result<Vec<Vec<f64>>, Error> {
    let n_atoms = masses.len();
    let size = 3 * n_atoms;
    let lattice = *structure.lattice();
    let cart: Vec<[f64; 3]> = structure
        .sites()
        .iter()
        .map(|site| lattice.to_cart(site.frac()))
        .collect();
    let total_mass: f64 = masses.iter().sum();
    if total_mass <= 0.0 {
        return Err(Error::InvalidInput(
            "a molecule or crystal needs a positive total mass".to_owned(),
        ));
    }
    let mut centre = [0.0; 3];
    for (mass, point) in masses.iter().zip(&cart) {
        centre[0] += mass * point[0];
        centre[1] += mass * point[1];
        centre[2] += mass * point[2];
    }
    for value in &mut centre {
        *value /= total_mass;
    }

    let mut candidates: Vec<Vec<f64>> = Vec::new();
    for axis in 0..3 {
        let mut vector = vec![0.0; size];
        for (index, mass) in masses.iter().enumerate() {
            if let Some(slot) = vector.get_mut(3 * index + axis) {
                *slot = mass.sqrt();
            }
        }
        candidates.push(vector);
    }
    if mode == DisplacementMode::Molecule {
        for axis in 0..3 {
            let mut unit = [0.0; 3];
            if let Some(slot) = unit.get_mut(axis) {
                *slot = 1.0;
            }
            let mut vector = vec![0.0; size];
            for (index, mass) in masses.iter().enumerate() {
                let Some(point) = cart.get(index) else {
                    continue;
                };
                let relative = [
                    point[0] - centre[0],
                    point[1] - centre[1],
                    point[2] - centre[2],
                ];
                let rotated = cross(unit, relative);
                for component in 0..3 {
                    if let Some(value) = rotated.get(component)
                        && let Some(slot) = vector.get_mut(3 * index + component)
                    {
                        *slot = mass.sqrt() * value;
                    }
                }
            }
            candidates.push(vector);
        }
    }

    let mut basis: Vec<Vec<f64>> = Vec::new();
    for candidate in candidates {
        let mut vector = candidate;
        for previous in &basis {
            let projection = dot(&vector, previous);
            for (value, base) in vector.iter_mut().zip(previous) {
                *value -= projection * base;
            }
        }
        let length = dot(&vector, &vector).sqrt();
        if length > 1.0e-10 {
            for value in &mut vector {
                *value /= length;
            }
            basis.push(vector);
        }
    }
    Ok(basis)
}

/// Projects `matrix` onto the orthogonal complement of `basis` on both sides.
pub(crate) fn project_out(matrix: &mut [Vec<f64>], basis: &[Vec<f64>]) {
    let n = matrix.len();
    for vector in basis {
        let mut product = vec![0.0; n];
        for row in 0..n {
            let mut total = 0.0;
            if let Some(values) = matrix.get(row) {
                for (column, value) in values.iter().enumerate() {
                    total += value * vector.get(column).copied().unwrap_or(0.0);
                }
            }
            if let Some(slot) = product.get_mut(row) {
                *slot = total;
            }
        }
        let rayleigh: f64 = vector
            .iter()
            .zip(&product)
            .map(|(first, second)| first * second)
            .sum();
        for row in 0..n {
            let first = vector.get(row).copied().unwrap_or(0.0);
            let second = product.get(row).copied().unwrap_or(0.0);
            if let Some(values) = matrix.get_mut(row) {
                for (column, value) in values.iter_mut().enumerate() {
                    let base = vector.get(column).copied().unwrap_or(0.0);
                    let other = product.get(column).copied().unwrap_or(0.0);
                    *value -= first * other + second * base - rayleigh * first * base;
                }
            }
        }
    }
}

fn dot(first: &[f64], second: &[f64]) -> f64 {
    first.iter().zip(second).map(|(a, b)| a * b).sum()
}
