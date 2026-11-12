#![allow(clippy::indexing_slicing)] // fixed 3x3 integer lattices

//! Integer lattice mappings for the structure matcher.

use xtalforge_core::{Error, Lattice};

use crate::MatcherOptions;

/// Finds an integer matrix `M` with `M * second ~ first` and determinant
/// `target_determinant`.
pub(crate) fn find_integer_mappings(
    first: &Lattice,
    second: &Lattice,
    options: &MatcherOptions,
    target_determinant: f64,
) -> Result<Vec<[[f64; 3]; 3]>, Error> {
    let inverse = second.inverse();
    let matrix = first.matrix();
    let mut candidates: Vec<Vec<[i32; 3]>> = Vec::new();
    for row in matrix {
        let fractional = [
            row[0] * inverse[0][0] + row[1] * inverse[1][0] + row[2] * inverse[2][0],
            row[0] * inverse[0][1] + row[1] * inverse[1][1] + row[2] * inverse[2][1],
            row[0] * inverse[0][2] + row[1] * inverse[1][2] + row[2] * inverse[2][2],
        ];
        let rounded = fractional.map(|value| value.round() as i32);
        let mut options_row = Vec::new();
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let candidate = [rounded[0] + dx, rounded[1] + dy, rounded[2] + dz];
                    let error = fractional
                        .iter()
                        .zip(candidate)
                        .map(|(value, integer)| (value - f64::from(integer)).abs())
                        .fold(0.0_f64, f64::max);
                    if error <= options.ltol {
                        options_row.push((error, candidate));
                    }
                }
            }
        }
        options_row.sort_by(|a, b| {
            a.0.partial_cmp(&b.0)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.1.cmp(&b.1))
        });
        options_row.dedup_by_key(|(_, candidate)| *candidate);
        if options_row.is_empty() {
            return Ok(Vec::new());
        }
        candidates.push(options_row.into_iter().map(|(_, vector)| vector).collect());
    }
    let mut mapped: Vec<([[i32; 3]; 3], f64)> = Vec::new();
    for a in candidates.first().map(Vec::as_slice).unwrap_or(&[]) {
        for b in candidates.get(1).map(Vec::as_slice).unwrap_or(&[]) {
            for c in candidates.get(2).map(Vec::as_slice).unwrap_or(&[]) {
                let matrix = [*a, *b, *c];
                let determinant = determinant_i32(matrix);
                let valid = if (target_determinant.abs() - 1.0).abs() < 1e-9 {
                    determinant.abs() == 1
                } else {
                    (determinant as f64 - target_determinant).abs() < 1e-6
                };
                if !valid {
                    continue;
                }
                let score = lattice_error(first, second, &matrix, options);
                if score <= options.ltol * 3.0 {
                    mapped.push((matrix, score));
                }
            }
        }
    }
    mapped.sort_by(|first_entry, second_entry| {
        first_entry
            .1
            .partial_cmp(&second_entry.1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| first_entry.0.cmp(&second_entry.0))
    });
    mapped.dedup_by_key(|(matrix, _)| *matrix);
    mapped.truncate(32);
    Ok(mapped
        .into_iter()
        .map(|(matrix, _)| matrix.map(|row| row.map(f64::from)))
        .collect())
}

fn lattice_error(
    first: &Lattice,
    second: &Lattice,
    matrix: &[[i32; 3]; 3],
    options: &MatcherOptions,
) -> f64 {
    let base = second.matrix();
    let lengths = first.lengths();
    let angles = first.angles();
    let mut error = 0.0;
    let mut transformed = [[0.0; 3]; 3];
    for row in 0..3 {
        for column in 0..3 {
            let mut sum = 0.0;
            for k in 0..3 {
                sum += f64::from(matrix[row][k]) * base[k][column];
            }
            transformed[row][column] = sum;
        }
    }
    for row in 0..3 {
        let length = (transformed[row][0] * transformed[row][0]
            + transformed[row][1] * transformed[row][1]
            + transformed[row][2] * transformed[row][2])
            .sqrt();
        if lengths[row] > 0.0 {
            error += ((length - lengths[row]) / lengths[row]).abs();
        }
    }
    for (index, (i, j)) in [(0, 1), (0, 2), (1, 2)].iter().enumerate() {
        let a = transformed[*i];
        let b = transformed[*j];
        let norm_a = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
        let norm_b = (b[0] * b[0] + b[1] * b[1] + b[2] * b[2]).sqrt();
        if norm_a < 1e-12 || norm_b < 1e-12 {
            continue;
        }
        let cosine =
            ((a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (norm_a * norm_b)).clamp(-1.0, 1.0);
        let angle = cosine.acos().to_degrees();
        if let Some(target) = angles.get(index) {
            error += ((angle - target) / options.angle_tol.max(1e-6)).abs() * 0.01;
        }
    }
    error
}

fn determinant_i32(matrix: [[i32; 3]; 3]) -> i64 {
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    i64::from(a) * (i64::from(e) * i64::from(i) - i64::from(f) * i64::from(h))
        - i64::from(b) * (i64::from(d) * i64::from(i) - i64::from(f) * i64::from(g))
        + i64::from(c) * (i64::from(d) * i64::from(h) - i64::from(e) * i64::from(g))
}
