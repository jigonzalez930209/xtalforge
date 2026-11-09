//! Wyckoff-allowed displacement directions for a site (Task 7.8.1).
//!
//! A site on a special position may only move along the directions its
//! stabilizer leaves invariant. For an operation `R` that maps the site onto
//! itself modulo a lattice translation, a displacement `u` keeps the site on
//! the same Wyckoff position when `R u = u`. Intersecting that condition over
//! every stabilizer element gives the free subspace: three directions on a
//! general position, one on a mirror or a rotation axis, none at an inversion
//! centre.

#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed 3x3 systems

use xtalforge_core::SymOp;

/// Tolerance for "this operation fixes the site".
const SITE_TOLERANCE: f64 = 1e-3;

/// Orthonormal basis of the allowed displacement directions.
///
/// The basis is empty at a fully fixed site and has up to three entries at a
/// general position, in a deterministic order.
#[must_use]
pub(crate) fn allowed_directions(operations: &[SymOp], frac: [f64; 3]) -> Vec<[f64; 3]> {
    // A = sum (R - I)^T (R - I) over the stabilizer; its null space is the
    // allowed subspace.
    let mut accumulator = [[0.0; 3]; 3];
    for operation in operations {
        if !same_point(operation.apply(frac), frac) {
            continue;
        }
        let rotation = operation.rotation();
        let mut delta = [[0.0; 3]; 3];
        for (row, line) in delta.iter_mut().enumerate() {
            for (column, value) in line.iter_mut().enumerate() {
                *value = f64::from(rotation[row][column]) - if row == column { 1.0 } else { 0.0 };
            }
        }
        for first in 0..3 {
            for second in 0..3 {
                let mut sum = 0.0;
                for middle in 0..3 {
                    sum += delta[middle][first] * delta[middle][second];
                }
                accumulator[first][second] += sum;
            }
        }
    }
    null_space(accumulator)
}

fn same_point(first: [f64; 3], second: [f64; 3]) -> bool {
    (0..3).all(|axis| {
        let delta =
            first.get(axis).copied().unwrap_or(0.0) - second.get(axis).copied().unwrap_or(0.0);
        let wrapped = delta - delta.round();
        wrapped.abs() <= SITE_TOLERANCE
    })
}

/// Orthonormal basis of the null space of a symmetric 3x3 matrix.
fn null_space(matrix: [[f64; 3]; 3]) -> Vec<[f64; 3]> {
    let scale = matrix
        .iter()
        .flat_map(|row| row.iter())
        .fold(0.0_f64, |value, entry| value.max(entry.abs()));
    if scale < 1e-12 {
        return vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }
    // Reduced row echelon form to find the free columns.
    let mut rows = matrix;
    let tolerance = 1e-9 * scale;
    let mut pivots = Vec::new();
    let mut row = 0;
    for column in 0..3 {
        let Some(pivot) = (row..3).max_by(|a, b| {
            rows[*a][column]
                .abs()
                .partial_cmp(&rows[*b][column].abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        }) else {
            continue;
        };
        if rows[pivot][column].abs() <= tolerance {
            continue;
        }
        rows.swap(row, pivot);
        let divisor = rows[row][column];
        for value in rows[row].iter_mut() {
            *value /= divisor;
        }
        for other in 0..3 {
            if other == row {
                continue;
            }
            let factor = rows[other][column];
            if factor == 0.0 {
                continue;
            }
            for index in 0..3 {
                rows[other][index] -= factor * rows[row][index];
            }
        }
        pivots.push(column);
        row += 1;
    }
    let mut directions = Vec::new();
    for free in 0..3 {
        if pivots.contains(&free) {
            continue;
        }
        let mut vector = [0.0; 3];
        if let Some(slot) = vector.get_mut(free) {
            *slot = 1.0;
        }
        for (index, pivot) in pivots.iter().enumerate() {
            if let Some(value) = rows.get(index).and_then(|line| line.get(free))
                && let Some(slot) = vector.get_mut(*pivot)
            {
                *slot = -*value;
            }
        }
        directions.push(normalize(vector));
    }
    orthonormalize(directions)
}

fn normalize(vector: [f64; 3]) -> [f64; 3] {
    let norm = vector.iter().map(|value| value * value).sum::<f64>().sqrt();
    if norm <= 1e-12 {
        [0.0; 3]
    } else {
        vector.map(|value| value / norm)
    }
}

fn orthonormalize(directions: Vec<[f64; 3]>) -> Vec<[f64; 3]> {
    let mut basis: Vec<[f64; 3]> = Vec::new();
    for direction in directions {
        let mut vector = direction;
        for existing in &basis {
            let projection = vector
                .iter()
                .zip(existing.iter())
                .map(|(a, b)| a * b)
                .sum::<f64>();
            for (slot, value) in vector.iter_mut().zip(existing.iter()) {
                *slot -= projection * value;
            }
        }
        let normalised = normalize(vector);
        if normalised.iter().any(|value| *value != 0.0) {
            basis.push(normalised);
        }
    }
    basis
}
