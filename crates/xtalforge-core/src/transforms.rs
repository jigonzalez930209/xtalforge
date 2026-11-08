//! Point-set transforms, finite cells and displacement ellipsoids (Task 1.9).

#![allow(clippy::indexing_slicing)] // small fixed-size matrix and vector math

use crate::{Error, Site, Structure, StructureBuilder};

/// A Cartesian plane through a point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Plane {
    /// Unit normal.
    pub normal: [f64; 3],
    /// A point on the plane.
    pub point: [f64; 3],
}

/// Displacement ellipsoid axes and semi-axis lengths.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ellipsoid {
    /// Three orthonormal axes.
    pub axes: [[f64; 3]; 3],
    /// Semi-axis lengths in Angstrom.
    pub lengths: [f64; 3],
    /// Probability level used for the lengths.
    pub probability: f64,
}

/// Translates Cartesian points.
#[must_use]
pub fn translate(points: &[[f64; 3]], delta: [f64; 3]) -> Vec<[f64; 3]> {
    points
        .iter()
        .map(|point| {
            [
                point[0] + delta[0],
                point[1] + delta[1],
                point[2] + delta[2],
            ]
        })
        .collect()
}

/// Rotates Cartesian points about an axis through `center`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a zero axis or non-finite values.
pub fn rotate(
    points: &[[f64; 3]],
    axis: [f64; 3],
    angle_degrees: f64,
    center: [f64; 3],
) -> Result<Vec<[f64; 3]>, Error> {
    if !axis.iter().all(|value| value.is_finite()) || !angle_degrees.is_finite() {
        return Err(Error::InvalidInput(
            "rotation axis and angle must be finite".to_owned(),
        ));
    }
    let axis = normalize(axis)?;
    let angle = angle_degrees.to_radians();
    let (sine, cosine) = angle.sin_cos();
    Ok(points
        .iter()
        .map(|point| {
            let relative = subtract(*point, center);
            let rotated = add(
                add(scale(relative, cosine), scale(cross(axis, relative), sine)),
                scale(axis, dot(axis, relative) * (1.0 - cosine)),
            );
            add(rotated, center)
        })
        .collect())
}

/// Mirrors Cartesian points through a plane.
#[must_use]
pub fn mirror(points: &[[f64; 3]], plane: Plane) -> Vec<[f64; 3]> {
    points
        .iter()
        .map(|point| {
            let distance = dot(subtract(*point, plane.point), plane.normal);
            subtract(*point, scale(plane.normal, 2.0 * distance))
        })
        .collect()
}

/// Inverts Cartesian points through a centre.
#[must_use]
pub fn invert(points: &[[f64; 3]], center: [f64; 3]) -> Vec<[f64; 3]> {
    points
        .iter()
        .map(|point| {
            [
                2.0 * center[0] - point[0],
                2.0 * center[1] - point[1],
                2.0 * center[2] - point[2],
            ]
        })
        .collect()
}

/// Rotates points so `from` aligns with `to`, about `center`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for zero or non-finite vectors.
pub fn align_vector(
    points: &[[f64; 3]],
    from: [f64; 3],
    to: [f64; 3],
    center: [f64; 3],
) -> Result<Vec<[f64; 3]>, Error> {
    let from = normalize(from)?;
    let to = normalize(to)?;
    let axis = cross(from, to);
    let sine = norm(axis);
    if sine < 1e-12 {
        if dot(from, to) > 0.0 {
            return Ok(points.to_vec());
        }
        let perpendicular = if from[0].abs() < 0.9 {
            normalize(cross(from, [1.0, 0.0, 0.0]))?
        } else {
            normalize(cross(from, [0.0, 1.0, 0.0]))?
        };
        return rotate(points, perpendicular, 180.0, center);
    }
    let angle = sine.atan2(dot(from, to)).to_degrees();
    rotate(points, axis, angle, center)
}

/// Fits a plane through at least three non-collinear points.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for fewer than three points or collinear
/// points.
pub fn fit_plane(points: &[[f64; 3]]) -> Result<Plane, Error> {
    if points.len() < 3 {
        return Err(Error::InvalidInput(
            "a plane needs at least three points".to_owned(),
        ));
    }
    let centroid = centroid(points);
    let mut first = [0.0; 3];
    for point in points {
        let candidate = subtract(*point, centroid);
        if norm(candidate) > 1e-9 {
            first = candidate;
            break;
        }
    }
    let mut second = [0.0; 3];
    for point in points {
        let candidate = subtract(*point, centroid);
        let normal = cross(first, candidate);
        if norm(normal) > 1e-9 {
            second = candidate;
            break;
        }
    }
    let normal = normalize(cross(first, second))?;
    Ok(Plane {
        normal,
        point: centroid,
    })
}

/// Rotates points so one fitted plane aligns with another.
///
/// # Errors
///
/// Returns the errors of [`fit_plane`] or [`align_vector`].
pub fn align_planes(
    points: &[[f64; 3]],
    from_normal: [f64; 3],
    to_normal: [f64; 3],
    center: [f64; 3],
) -> Result<Vec<[f64; 3]>, Error> {
    align_vector(points, from_normal, to_normal, center)
}

/// Converts anisotropic displacement parameters to an ellipsoid.
///
/// `uij` is `[U11, U22, U33, U12, U13, U23]` in Angstrom squared. The
/// probability is one of 0.5, 0.9, 0.95 or 0.99; other values are interpolated
/// from the chi-square table for three degrees of freedom.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for non-finite or negative-definite tensors.
pub fn displacement_ellipsoid(uij: [f64; 6], probability: f64) -> Result<Ellipsoid, Error> {
    if !uij.iter().all(|value| value.is_finite()) || !probability.is_finite() {
        return Err(Error::InvalidInput(
            "displacement parameters must be finite".to_owned(),
        ));
    }
    if !(0.0 < probability && probability < 1.0) {
        return Err(Error::InvalidInput(
            "ellipsoid probability must be in (0, 1)".to_owned(),
        ));
    }
    let matrix = [
        [uij[0], uij[3], uij[4]],
        [uij[3], uij[1], uij[5]],
        [uij[4], uij[5], uij[2]],
    ];
    let (values, vectors) = jacobi_eigen(matrix)?;
    let mut lengths = [0.0; 3];
    for index in 0..3 {
        if values[index] < -1e-12 {
            return Err(Error::InvalidInput(
                "displacement tensor is not positive semidefinite".to_owned(),
            ));
        }
        lengths[index] = values[index].max(0.0).sqrt() * chi2_scale(probability);
    }
    Ok(Ellipsoid {
        axes: vectors,
        lengths,
        probability,
    })
}

impl Structure {
    /// Atoms of one finite cell, including images on faces, edges and corners.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a negative tolerance.
    pub fn finite_cell(&self, tolerance: f64) -> Result<Self, Error> {
        if !tolerance.is_finite() || tolerance < 0.0 {
            return Err(Error::InvalidInput(
                "finite-cell tolerance must be non-negative and finite".to_owned(),
            ));
        }
        let p1 = self.p1()?;
        if p1.atoms.is_empty() {
            return Ok(self.clone());
        }
        let lattice = p1.lattice;
        let mut sites: Vec<Site> = Vec::new();
        for atom in &p1.atoms {
            let shifts = atom.frac.map(|value| {
                let mut out = vec![0];
                if value <= tolerance {
                    out.push(1);
                }
                if value >= 1.0 - tolerance {
                    out.push(-1);
                }
                out
            });
            for x in &shifts[0] {
                for y in &shifts[1] {
                    for z in &shifts[2] {
                        let frac = [
                            atom.frac[0] + f64::from(*x),
                            atom.frac[1] + f64::from(*y),
                            atom.frac[2] + f64::from(*z),
                        ];
                        if frac
                            .iter()
                            .all(|value| (-tolerance..=1.0 + tolerance).contains(value))
                        {
                            let duplicate = sites.iter().any(|existing| {
                                existing.element() == atom.element
                                    && existing.frac().iter().zip(frac.iter()).all(
                                        |(first, second)| {
                                            (first - second).abs() <= tolerance.max(1e-9)
                                        },
                                    )
                            });
                            if !duplicate {
                                sites.push(Site::new(
                                    format!("{}{}", atom.element, sites.len() + 1),
                                    &atom.element,
                                    frac,
                                )?);
                            }
                        }
                    }
                }
            }
        }
        let mut builder = StructureBuilder::new(self.name(), lattice);
        for site in sites {
            builder = builder.site(site);
        }
        builder.build()
    }
}

fn centroid(points: &[[f64; 3]]) -> [f64; 3] {
    let mut sum = [0.0; 3];
    for point in points {
        sum[0] += point[0];
        sum[1] += point[1];
        sum[2] += point[2];
    }
    let count = points.len() as f64;
    [sum[0] / count, sum[1] / count, sum[2] / count]
}

fn chi2_scale(probability: f64) -> f64 {
    let table = [(0.5, 1.538), (0.9, 2.500), (0.95, 2.796), (0.99, 3.368)];
    if probability <= table[0].0 {
        return table[0].1;
    }
    let mut previous = table[0];
    for (level, value) in table {
        if probability <= level {
            let fraction = (probability - previous.0) / (level - previous.0);
            return previous.1 + fraction * (value - previous.1);
        }
        previous = (level, value);
    }
    table[table.len() - 1].1
}

fn jacobi_eigen(matrix: [[f64; 3]; 3]) -> Result<([f64; 3], [[f64; 3]; 3]), Error> {
    let mut a = matrix;
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..32 {
        let (p, q) = largest_off_diagonal(a);
        if a[p][q].abs() < 1e-12 {
            break;
        }
        let angle = 0.5 * (2.0 * a[p][q]).atan2(a[q][q] - a[p][p]);
        let (sine, cosine) = angle.sin_cos();
        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];
        a[p][p] = cosine * cosine * app - 2.0 * sine * cosine * apq + sine * sine * aqq;
        a[q][q] = sine * sine * app + 2.0 * sine * cosine * apq + cosine * cosine * aqq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;
        for index in 0..3 {
            if index != p && index != q {
                let aip = a[index][p];
                let aiq = a[index][q];
                a[index][p] = cosine * aip - sine * aiq;
                a[p][index] = a[index][p];
                a[index][q] = sine * aip + cosine * aiq;
                a[q][index] = a[index][q];
            }
            let vpi = v[index][p];
            let vqi = v[index][q];
            v[index][p] = cosine * vpi - sine * vqi;
            v[index][q] = sine * vpi + cosine * vqi;
        }
    }
    let mut values = [a[0][0], a[1][1], a[2][2]];
    let mut axes = [
        [v[0][0], v[1][0], v[2][0]],
        [v[0][1], v[1][1], v[2][1]],
        [v[0][2], v[1][2], v[2][2]],
    ];
    for index in 0..3 {
        if values[index] < -1e-9 {
            return Err(Error::InvalidInput(
                "displacement tensor is not positive semidefinite".to_owned(),
            ));
        }
        values[index] = values[index].max(0.0);
        axes[index] = normalize(axes[index]).unwrap_or([0.0, 0.0, 0.0]);
    }
    let mut order = [0usize, 1, 2];
    order.sort_by(|first, second| {
        values[*first]
            .partial_cmp(&values[*second])
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok((
        [values[order[0]], values[order[1]], values[order[2]]],
        [axes[order[0]], axes[order[1]], axes[order[2]]],
    ))
}

fn largest_off_diagonal(matrix: [[f64; 3]; 3]) -> (usize, usize) {
    let candidates = [(0, 1), (0, 2), (1, 2)];
    let mut best = candidates[0];
    for candidate in candidates {
        if matrix[candidate.0][candidate.1].abs() > matrix[best.0][best.1].abs() {
            best = candidate;
        }
    }
    best
}

fn add(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[0] + second[0],
        first[1] + second[1],
        first[2] + second[2],
    ]
}

fn subtract(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[0] - second[0],
        first[1] - second[1],
        first[2] - second[2],
    ]
}

fn scale(vector: [f64; 3], factor: f64) -> [f64; 3] {
    [vector[0] * factor, vector[1] * factor, vector[2] * factor]
}

fn cross(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    ]
}

fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    first[0] * second[0] + first[1] * second[1] + first[2] * second[2]
}

fn norm(vector: [f64; 3]) -> f64 {
    dot(vector, vector).sqrt()
}

fn normalize(vector: [f64; 3]) -> Result<[f64; 3], Error> {
    let length = norm(vector);
    if length < 1e-12 {
        return Err(Error::InvalidInput("the vector has zero length".to_owned()));
    }
    Ok(scale(vector, 1.0 / length))
}
