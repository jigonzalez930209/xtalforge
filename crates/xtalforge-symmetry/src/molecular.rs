//! Molecular point groups for finite molecules (Task 2.8.3).
//!
//! The detector builds candidate rotations by mapping a non-collinear triple
//! of atom vectors onto every triple with the same element pattern and Gram
//! matrix, then keeps the candidates that permute the whole molecule. The
//! Schoenflies symbol is read from the counts of proper rotations, mirrors and
//! inversion. Linear molecules (infinite groups) are rejected.

use nalgebra::{Matrix3, Vector3};
use serde::Serialize;
use xtalforge_core::Error;

/// Point-group summary of a molecule.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct MolecularPointGroup {
    /// Schoenflies symbol, for example `C2v` or `Td`.
    pub schoenflies: String,
    /// Group order, including improper operations.
    pub order: usize,
    /// Number of proper rotations.
    pub proper_rotations: usize,
    /// Number of improper operations.
    pub improper_operations: usize,
    /// Whether the inversion is present.
    pub has_inversion: bool,
    /// Number of mirror planes.
    pub mirrors: usize,
}

/// Detects the point group of a finite molecule.
///
/// The points are Cartesian coordinates in Angstrom; elements are used to
/// restrict candidate permutations.
///
/// # Errors
///
/// Returns [`Error::Symmetry`] for fewer than two atoms, a collinear molecule
/// or a group outside the supported Schoenflies table.
pub fn molecular_point_group(
    points: &[(String, [f64; 3])],
    tolerance: f64,
) -> Result<MolecularPointGroup, Error> {
    if points.len() < 2 {
        return Err(Error::Symmetry(
            "a point group needs at least two atoms".to_owned(),
        ));
    }
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(Error::Symmetry(format!(
            "point-group tolerance must be positive, got {tolerance}"
        )));
    }
    let count = points.len() as f64;
    let mut center = Vector3::zeros();
    for (_, point) in points {
        center += Vector3::new(point[0], point[1], point[2]);
    }
    center /= count;
    let vectors: Vec<Vector3<f64>> = points
        .iter()
        .map(|(_, point)| Vector3::new(point[0], point[1], point[2]) - center)
        .collect();
    let reference = reference_pair(&vectors, tolerance)?;
    let mut operations: Vec<Matrix3<f64>> = Vec::new();
    for first in 0..points.len() {
        for second in 0..points.len() {
            if first == second {
                continue;
            }
            let target = [first, second];
            if !same_pattern(points, &reference, &target) {
                continue;
            }
            for sign in [1.0, -1.0] {
                let Some(candidate) = build_rotation(&vectors, &reference, &target, sign) else {
                    continue;
                };
                if operations
                    .iter()
                    .any(|existing| (existing - candidate).norm() < 1e-4)
                {
                    continue;
                }
                if maps_molecule(&vectors, &candidate, points, tolerance) {
                    operations.push(candidate);
                }
            }
        }
    }
    classify(&operations)
}

fn reference_pair(vectors: &[Vector3<f64>], tolerance: f64) -> Result<[usize; 2], Error> {
    for first in 0..vectors.len() {
        for second in 0..vectors.len() {
            let (Some(first_vector), Some(second_vector)) =
                (vectors.get(first), vectors.get(second))
            else {
                continue;
            };
            if first != second && first_vector.cross(second_vector).norm() > tolerance {
                return Ok([first, second]);
            }
        }
    }
    Err(Error::Symmetry(
        "the molecule has no two non-collinear atoms".to_owned(),
    ))
}

fn same_pattern(
    points: &[(String, [f64; 3])],
    reference: &[usize; 2],
    target: &[usize; 2],
) -> bool {
    reference.iter().zip(target.iter()).all(|(first, second)| {
        points
            .get(*first)
            .zip(points.get(*second))
            .is_some_and(|(first, second)| first.0 == second.0)
    })
}

fn build_rotation(
    vectors: &[Vector3<f64>],
    reference: &[usize; 2],
    target: &[usize; 2],
    sign: f64,
) -> Option<Matrix3<f64>> {
    let source = frame(*vectors.get(reference[0])?, *vectors.get(reference[1])?)?;
    let destination = frame(*vectors.get(target[0])?, *vectors.get(target[1])?)?;
    let destination = destination * Matrix3::from_diagonal(&Vector3::new(1.0, 1.0, sign));
    let candidate = destination * source.try_inverse()?;
    let identity = Matrix3::identity();
    if (candidate.transpose() * candidate - identity).norm() > 1e-3 {
        return None;
    }
    Some(candidate)
}

fn frame(first: Vector3<f64>, second: Vector3<f64>) -> Option<Matrix3<f64>> {
    let first = first.try_normalize(1e-12)?;
    let normal = first.cross(&second).try_normalize(1e-12)?;
    let second = normal.cross(&first);
    Some(Matrix3::from_columns(&[first, second, normal]))
}

fn maps_molecule(
    vectors: &[Vector3<f64>],
    rotation: &Matrix3<f64>,
    points: &[(String, [f64; 3])],
    tolerance: f64,
) -> bool {
    for (index, (element, _)) in points.iter().enumerate() {
        let Some(vector) = vectors.get(index) else {
            return false;
        };
        let image = rotation * vector;
        let matched = points
            .iter()
            .enumerate()
            .any(|(other, (other_element, _))| {
                other_element == element
                    && vectors
                        .get(other)
                        .is_some_and(|candidate| (candidate - image).norm() <= tolerance)
            });
        if !matched {
            return false;
        }
    }
    true
}

fn classify(operations: &[Matrix3<f64>]) -> Result<MolecularPointGroup, Error> {
    let mut proper = 0;
    let mut improper = 0;
    let mut has_inversion = false;
    let mut mirrors = 0;
    let mut rotation_orders: Vec<usize> = Vec::new();
    for operation in operations {
        let determinant = operation.determinant();
        if determinant > 0.0 {
            proper += 1;
            let cosine = ((operation.trace() - 1.0) / 2.0).clamp(-1.0, 1.0);
            let angle = cosine.acos();
            if angle > 1e-6 {
                let order = (std::f64::consts::TAU / angle).round() as usize;
                rotation_orders.push(order.max(2));
            }
        } else {
            improper += 1;
            if operation.trace() < -2.5 {
                has_inversion = true;
            } else if operation.trace() > 0.5 && operation.trace() < 1.5 {
                mirrors += 1;
            }
        }
    }
    let max_order = rotation_orders.iter().copied().max().unwrap_or(1);
    let schoenflies = match (
        operations.len(),
        proper,
        improper,
        has_inversion,
        mirrors,
        max_order,
    ) {
        (1, _, _, _, _, _) => "C1",
        (2, 1, 1, _, 0, _) => "Ci",
        (2, 1, 1, _, 1, _) => "Cs",
        (2, 2, 0, _, _, 2) => "C2",
        (4, 2, 2, _, 2, 2) => "C2v",
        (4, 2, 2, true, 1, 2) => "C2h",
        (4, 4, 0, _, _, 2) => "D2",
        (8, 4, 4, true, 3, 2) => "D2h",
        (3, 3, 0, _, _, 3) => "C3",
        (6, 3, 3, _, 3, 3) => "C3v",
        (6, 3, 3, true, 0, 3) => "C3i",
        (6, 6, 0, _, _, 2) => "D3",
        (12, 6, 6, true, 3, 3) => "D3d",
        (12, 6, 6, false, 4, 3) => "D3h",
        (4, 4, 0, _, _, 4) => "C4",
        (8, 4, 4, _, 4, 4) => "C4v",
        (8, 4, 4, true, 1, 4) => "C4h",
        (8, 8, 0, _, _, 2) => "D4",
        (16, 8, 8, true, 5, 4) => "D4h",
        (4, 2, 2, _, 0, 2) => "S4",
        (5, 5, 0, _, _, 5) => "C5",
        (10, 5, 5, _, 5, 5) => "C5v",
        (10, 10, 0, _, _, 2) => "D5",
        (20, 10, 10, _, 6, 5) => "D5h",
        (6, 6, 0, _, _, 6) => "C6",
        (12, 6, 6, _, 6, 6) => "C6v",
        (12, 12, 0, _, _, 2) => "D6",
        (24, 12, 12, true, 7, 6) => "D6h",
        (12, 12, 0, _, _, 3) => "T",
        (24, 12, 12, _, 6, 3) => "Td",
        (24, 12, 12, true, 3, 3) => "Th",
        (24, 24, 0, _, _, 4) => "O",
        (48, 24, 24, true, 9, 4) => "Oh",
        _ => {
            return Err(Error::Symmetry(format!(
                "unsupported molecular point group: order {}, {} proper, {} improper, {} mirrors, max order {}",
                operations.len(),
                proper,
                improper,
                mirrors,
                max_order
            )));
        }
    };
    Ok(MolecularPointGroup {
        schoenflies: schoenflies.to_owned(),
        order: operations.len(),
        proper_rotations: proper,
        improper_operations: improper,
        has_inversion,
        mirrors,
    })
}
