//! Basis-independent fingerprint of an ordered structure.
//!
//! [`match_structures`](xtalforge_symmetry::match_structures) handles most
//! duplicates, but its primitive-cell reduction is not guaranteed for every
//! two-atom cell. This fallback canonicalizes the standardized cell under the
//! lattice point group, which keeps the derivative list free of duplicates
//! even when the matcher declines a pair.

use xtalforge_core::{Lattice, Structure};
use xtalforge_symmetry::{SymmetryDataset, SymmetryOptions, detect_symmetry};

/// Canonical fingerprint of one structure.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct CanonicalSignature {
    /// International space-group number.
    number: u32,
    /// Sorted, quantized metric-tensor entries.
    metric: Vec<i64>,
    /// Sorted `(element, x, y, z)` positions after the origin shift.
    sites: Vec<(String, i64, i64, i64)>,
}

/// Returns the canonical signature and the number of symmetry operations.
///
/// The signature is absent when the symmetry search fails; the operation count
/// is then one.
pub(super) fn signature_and_symmetry(structure: &Structure) -> (Option<CanonicalSignature>, usize) {
    let dataset = detect_symmetry(
        structure,
        SymmetryOptions {
            symprec: 0.05,
            ..SymmetryOptions::default()
        },
    );
    let Ok(dataset) = dataset else {
        return (None, 1);
    };
    let count = dataset.operations.len();
    let signature = signature_from_dataset(&dataset);
    (signature, count)
}

fn signature_from_dataset(dataset: &SymmetryDataset) -> Option<CanonicalSignature> {
    let lattice = Lattice::new(dataset.std_lattice).ok()?;
    let metric = lattice.metric_tensor();
    let mut metric_key: Vec<i64> = metric
        .iter()
        .flatten()
        .map(|value| (value * 1_000_000.0).round() as i64)
        .collect();
    metric_key.sort_unstable();
    let sites = canonical_sites(&dataset.std_positions, &dataset.std_elements, &metric)?;
    Some(CanonicalSignature {
        number: dataset.number,
        metric: metric_key,
        sites,
    })
}

fn canonical_sites(
    positions: &[[f64; 3]],
    elements: &[String],
    metric: &[[f64; 3]; 3],
) -> Option<Vec<(String, i64, i64, i64)>> {
    let mut best: Option<Vec<(String, i64, i64, i64)>> = None;
    for rotation in lattice_automorphisms(metric) {
        let rotation = rotation.map(|row| row.map(f64::from));
        let mut atoms: Vec<(String, [f64; 3])> = elements
            .iter()
            .zip(positions)
            .map(|(element, position)| {
                (
                    element.clone(),
                    Lattice::wrap(multiply_row(*position, rotation)),
                )
            })
            .collect();
        if atoms.is_empty() {
            return None;
        }
        let mut pivot: Option<(usize, (String, [i64; 3]))> = None;
        for (index, (element, position)) in atoms.iter().enumerate() {
            let key = (element.clone(), position_key(*position));
            if pivot.as_ref().is_none_or(|(_, current)| &key < current) {
                pivot = Some((index, key));
            }
        }
        let Some((index, _)) = pivot else {
            continue;
        };
        let shift = atoms.get(index).map(|(_, position)| *position)?;
        let mut signature: Vec<(String, i64, i64, i64)> = atoms
            .drain(..)
            .map(|(element, position)| {
                let shifted = Lattice::wrap([
                    position[0] - shift[0],
                    position[1] - shift[1],
                    position[2] - shift[2],
                ]);
                let [x, y, z] = position_key(shifted);
                (element, x, y, z)
            })
            .collect();
        signature.sort();
        if best.as_ref().is_none_or(|current| &signature < current) {
            best = Some(signature);
        }
    }
    best
}

fn lattice_automorphisms(metric: &[[f64; 3]; 3]) -> Vec<[[i32; 3]; 3]> {
    let mut operations = Vec::new();
    for a in -1..=1 {
        for b in -1..=1 {
            for c in -1..=1 {
                for d in -1..=1 {
                    for e in -1..=1 {
                        for f in -1..=1 {
                            for g in -1..=1 {
                                for h in -1..=1 {
                                    for i in -1..=1 {
                                        let rotation = [[a, b, c], [d, e, f], [g, h, i]];
                                        if determinant(rotation).abs() != 1 {
                                            continue;
                                        }
                                        if metric_preserved(rotation, metric) {
                                            operations.push(rotation);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    operations
}

fn metric_preserved(rotation: [[i32; 3]; 3], metric: &[[f64; 3]; 3]) -> bool {
    let rotation = rotation.map(|row| row.map(f64::from));
    let product = multiply(rotation, *metric);
    let transformed = multiply(product, transpose(rotation));
    let scale = metric
        .iter()
        .flatten()
        .fold(0.0_f64, |current, value| current.max(value.abs()))
        .max(1.0);
    transformed
        .iter()
        .flatten()
        .zip(metric.iter().flatten())
        .all(|(left, right)| (left - right).abs() <= 1e-6 * scale)
}

fn determinant(matrix: [[i32; 3]; 3]) -> i32 {
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
}

fn multiply(left: [[f64; 3]; 3], right: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let [[a, b, c], [d, e, f], [g, h, i]] = left;
    let [[j, k, l], [m, n, o], [p, q, r]] = right;
    [
        [
            a * j + b * m + c * p,
            a * k + b * n + c * q,
            a * l + b * o + c * r,
        ],
        [
            d * j + e * m + f * p,
            d * k + e * n + f * q,
            d * l + e * o + f * r,
        ],
        [
            g * j + h * m + i * p,
            g * k + h * n + i * q,
            g * l + h * o + i * r,
        ],
    ]
}

fn transpose(matrix: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    [[a, d, g], [b, e, h], [c, f, i]]
}

fn multiply_row(vector: [f64; 3], matrix: [[f64; 3]; 3]) -> [f64; 3] {
    let [x, y, z] = vector;
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    [
        x * a + y * d + z * g,
        x * b + y * e + z * h,
        x * c + y * f + z * i,
    ]
}

fn position_key(position: [f64; 3]) -> [i64; 3] {
    let wrapped = Lattice::wrap(position);
    wrapped.map(|value| (value * 1_000_000.0).round() as i64)
}
