#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Comparison primitives: Kabsch alignment, periodic RMSD and interpolation.
//!
//! The structure matcher of Task 9.1 lives in the symmetry crate, where the
//! Niggli reduction is available. This module holds the geometry that the
//! matcher, the binding surfaces and the interpolation share.

use serde::{Deserialize, Serialize};

use crate::{Error, Lattice, Structure};

/// A rigid alignment of one point set onto another.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Alignment {
    /// Rotation matrix in Cartesian coordinates.
    pub rotation: [[f64; 3]; 3],
    /// Translation applied after the rotation.
    pub translation: [f64; 3],
    /// Root mean square displacement in Angstrom.
    pub rmsd: f64,
    /// The first point set after the rotation and translation.
    pub aligned: Vec<[f64; 3]>,
}

/// Optimal rigid alignment of `points` onto `reference` (Kabsch, 1976).
///
/// The algorithm uses Horn's quaternion form, so no external linear-algebra
/// routine is needed and the result is deterministic.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for fewer than two points, different
/// lengths or non-finite coordinates.
pub fn kabsch_alignment(points: &[[f64; 3]], reference: &[[f64; 3]]) -> Result<Alignment, Error> {
    if points.len() != reference.len() {
        return Err(Error::InvalidInput(format!(
            "alignment needs the same number of points, got {} and {}",
            points.len(),
            reference.len()
        )));
    }
    if points.len() < 2 {
        return Err(Error::InvalidInput(
            "alignment needs at least two points".to_owned(),
        ));
    }
    if !points
        .iter()
        .chain(reference)
        .flatten()
        .all(|value| value.is_finite())
    {
        return Err(Error::InvalidInput(
            "alignment points must be finite".to_owned(),
        ));
    }
    let first_centre = centroid(points);
    let second_centre = centroid(reference);
    let mut covariance = [[0.0; 3]; 3];
    for (point, target) in points.iter().zip(reference) {
        let a = subtract(*point, first_centre);
        let b = subtract(*target, second_centre);
        for row in 0..3 {
            for column in 0..3 {
                if let (Some(value), Some(other)) = (a.get(row), b.get(column))
                    && let Some(entry) = covariance.get_mut(row).and_then(|r| r.get_mut(column))
                {
                    *entry += value * other;
                }
            }
        }
    }
    let rotation = rotation_from_covariance(covariance);
    let rotated_centre = apply(rotation, first_centre);
    let translation = [
        second_centre[0] - rotated_centre[0],
        second_centre[1] - rotated_centre[1],
        second_centre[2] - rotated_centre[2],
    ];
    let mut aligned = Vec::with_capacity(points.len());
    let mut squared = 0.0;
    for (point, target) in points.iter().zip(reference) {
        let moved = apply(rotation, *point);
        let moved = [
            moved[0] + translation[0],
            moved[1] + translation[1],
            moved[2] + translation[2],
        ];
        let delta = subtract(moved, *target);
        squared += dot(delta, delta);
        aligned.push(moved);
    }
    Ok(Alignment {
        rotation,
        translation,
        rmsd: (squared / points.len() as f64).sqrt(),
        aligned,
    })
}

/// Periodic RMS displacement of two matched P1 cells.
///
/// `mapping[i]` is the index in `second` that corresponds to P1 atom `i` of
/// `first`. Each pair is measured with the minimum image of the first cell.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an incomplete mapping, different
/// lattice volumes or an empty mapping.
pub fn periodic_rmsd(
    first: &Structure,
    second: &Structure,
    mapping: &[usize],
) -> Result<f64, Error> {
    if mapping.is_empty() {
        return Err(Error::InvalidInput(
            "periodic RMSD needs at least one mapped atom".to_owned(),
        ));
    }
    if !first.lattice().almost_equal(second.lattice(), 1e-6) {
        return Err(Error::InvalidInput(
            "periodic RMSD needs the same cell; use the structure matcher first".to_owned(),
        ));
    }
    let first_cell = first.p1()?;
    let second_cell = second.p1()?;
    if mapping.len() != first_cell.atoms.len() {
        return Err(Error::InvalidInput(format!(
            "mapping has {} entries for {} atoms",
            mapping.len(),
            first_cell.atoms.len()
        )));
    }
    let mut squared = 0.0;
    for (atom, &target) in first_cell.atoms.iter().zip(mapping) {
        let Some(other) = second_cell.atoms.get(target) else {
            return Err(Error::InvalidInput(format!(
                "mapping index {target} is out of range"
            )));
        };
        let delta = subtract(other.frac, atom.frac);
        let (_, distance) = first.lattice().minimum_image(delta);
        squared += distance * distance;
    }
    Ok((squared / mapping.len() as f64).sqrt())
}

/// Linear interpolation of two structures with the same atom order.
///
/// `images` intermediate frames are returned; the endpoints are not included.
/// Fractional differences take the minimum image across the cell boundary.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for `images == 0`, different atom counts or
/// elements, or an inconsistent cell.
pub fn interpolate(
    first: &Structure,
    second: &Structure,
    images: u32,
) -> Result<Vec<Structure>, Error> {
    if images == 0 {
        return Err(Error::InvalidInput(
            "interpolation needs at least one image".to_owned(),
        ));
    }
    let start = first.p1()?;
    let end = second.p1()?;
    if start.atoms.len() != end.atoms.len() {
        return Err(Error::InvalidInput(format!(
            "interpolation needs the same atom count, got {} and {}",
            start.atoms.len(),
            end.atoms.len()
        )));
    }
    for (a, b) in start.atoms.iter().zip(&end.atoms) {
        if a.element != b.element {
            return Err(Error::InvalidInput(format!(
                "atom {} is {} in the first structure and {} in the second",
                a.label, a.element, b.element
            )));
        }
    }
    let matrix = start.lattice.matrix();
    let end_matrix = end.lattice.matrix();
    let mut frames = Vec::with_capacity(images as usize);
    for frame in 1..=images {
        let t = f64::from(frame) / f64::from(images + 1);
        let mut lattice_rows = [[0.0; 3]; 3];
        for row in 0..3 {
            for column in 0..3 {
                let a = matrix[row][column];
                let b = end_matrix[row][column];
                lattice_rows[row][column] = a + (b - a) * t;
            }
        }
        let lattice = Lattice::new(lattice_rows)?;
        let mut builder = crate::StructureBuilder::new(first.name(), lattice);
        for (a, b) in start.atoms.iter().zip(&end.atoms) {
            let (delta, _) = start.lattice.minimum_image(subtract(b.frac, a.frac));
            let frac = [
                a.frac[0] + delta[0] * t,
                a.frac[1] + delta[1] * t,
                a.frac[2] + delta[2] * t,
            ];
            let site = crate::Site::new(a.label.clone(), &a.element, Lattice::wrap(frac))?;
            builder = builder.site(site.with_occupancy(a.occupancy)?);
        }
        frames.push(builder.build()?);
    }
    Ok(frames)
}

/// Rotation matrix from Horn's 4x4 quaternion matrix.
fn rotation_from_covariance(s: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let [[sxx, sxy, sxz], [syx, syy, syz], [szx, szy, szz]] = s;
    let matrix = [
        [sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
        [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
        [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
        [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
    ];
    let (values, vectors) = jacobi4(matrix);
    let mut best = 0;
    for index in 1..4 {
        if values.get(index).copied().unwrap_or(f64::NEG_INFINITY)
            > values.get(best).copied().unwrap_or(f64::NEG_INFINITY)
        {
            best = index;
        }
    }
    let quaternion = [
        vectors[0][best],
        vectors[1][best],
        vectors[2][best],
        vectors[3][best],
    ];
    quaternion_matrix(quaternion)
}

fn jacobi4(mut matrix: [[f64; 4]; 4]) -> ([f64; 4], [[f64; 4]; 4]) {
    let mut vectors = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];
    for _ in 0..64 {
        let mut p = 0;
        let mut q = 1;
        let mut largest = 0.0;
        for row in 0..4 {
            for column in (row + 1)..4 {
                let value = matrix[row][column].abs();
                if value > largest {
                    largest = value;
                    p = row;
                    q = column;
                }
            }
        }
        if largest < 1e-14 {
            break;
        }
        let theta = (matrix[q][q] - matrix[p][p]) / (2.0 * matrix[p][q]);
        let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
        let c = 1.0 / (t * t + 1.0).sqrt();
        let s = t * c;
        let mut next = matrix;
        for k in 0..4 {
            next[p][k] = c * matrix[p][k] - s * matrix[q][k];
            next[q][k] = s * matrix[p][k] + c * matrix[q][k];
        }
        matrix = next;
        let mut transposed = matrix;
        for k in 0..4 {
            transposed[k][p] = c * matrix[k][p] - s * matrix[k][q];
            transposed[k][q] = s * matrix[k][p] + c * matrix[k][q];
        }
        matrix = transposed;
        for k in 0..4 {
            let a = vectors[k][p];
            let b = vectors[k][q];
            vectors[k][p] = c * a - s * b;
            vectors[k][q] = s * a + c * b;
        }
    }
    (
        [matrix[0][0], matrix[1][1], matrix[2][2], matrix[3][3]],
        vectors,
    )
}

fn quaternion_matrix(quaternion: [f64; 4]) -> [[f64; 3]; 3] {
    let norm = (quaternion[0] * quaternion[0]
        + quaternion[1] * quaternion[1]
        + quaternion[2] * quaternion[2]
        + quaternion[3] * quaternion[3])
        .sqrt();
    if norm < 1e-12 {
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }
    let [w, x, y, z] = quaternion.map(|value| value / norm);
    [
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y - z * w),
            2.0 * (x * z + y * w),
        ],
        [
            2.0 * (x * y + z * w),
            1.0 - 2.0 * (x * x + z * z),
            2.0 * (y * z - x * w),
        ],
        [
            2.0 * (x * z - y * w),
            2.0 * (y * z + x * w),
            1.0 - 2.0 * (x * x + y * y),
        ],
    ]
}

fn centroid(points: &[[f64; 3]]) -> [f64; 3] {
    let mut total = [0.0; 3];
    for point in points {
        total[0] += point[0];
        total[1] += point[1];
        total[2] += point[2];
    }
    let count = points.len().max(1) as f64;
    [total[0] / count, total[1] / count, total[2] / count]
}

fn apply(matrix: [[f64; 3]; 3], point: [f64; 3]) -> [f64; 3] {
    [
        matrix[0][0] * point[0] + matrix[0][1] * point[1] + matrix[0][2] * point[2],
        matrix[1][0] * point[0] + matrix[1][1] * point[1] + matrix[1][2] * point[2],
        matrix[2][0] * point[0] + matrix[2][1] * point[1] + matrix[2][2] * point[2],
    ]
}

fn subtract(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[0] - second[0],
        first[1] - second[1],
        first[2] - second[2],
    ]
}

fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    first[0] * second[0] + first[1] * second[1] + first[2] * second[2]
}
