//! Ideal-geometry direction construction for hydrogen placement.

use std::f64::consts::{PI, TAU};

use xtalforge_core::Error;

use super::DEGENERATE;

/// Directions for `missing` hydrogens around `neighbours`.
pub(super) fn hydrogen_directions(
    neighbours: &[[f64; 3]],
    missing: usize,
    angle_degrees: f64,
) -> Result<Vec<[f64; 3]>, Error> {
    let theta = angle_degrees.to_radians();
    match neighbours {
        [] => no_neighbour_directions(missing),
        [first] => one_neighbour_directions(*first, missing, theta),
        [first, second] => two_neighbour_directions(*first, *second, missing, theta),
        [first, second, third] => {
            let sum = add(add(*first, *second), *third);
            if missing != 1 {
                return Err(Error::InvalidInput(format!(
                    "cannot place {missing} hydrogens on an atom with three neighbours"
                )));
            }
            Ok(vec![normalized(scale(sum, -1.0))?])
        }
        _ => Err(Error::InvalidInput(format!(
            "cannot place {missing} hydrogens on an atom with {} neighbours",
            neighbours.len()
        ))),
    }
}

/// Directions for an isolated atom: the cardinal axes, tetrahedral for four.
fn no_neighbour_directions(missing: usize) -> Result<Vec<[f64; 3]>, Error> {
    let axes = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    match missing {
        0 => Ok(Vec::new()),
        1..=3 => Ok(axes.iter().take(missing).copied().collect()),
        4 => {
            let tetrahedron = [
                [1.0, 1.0, 1.0],
                [1.0, -1.0, -1.0],
                [-1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0],
            ];
            tetrahedron
                .iter()
                .map(|direction| normalized(*direction))
                .collect()
        }
        _ => Err(Error::InvalidInput(format!(
            "cannot place {missing} hydrogens on an isolated atom"
        ))),
    }
}

/// Directions on a cone around the single existing bond.
fn one_neighbour_directions(
    first: [f64; 3],
    missing: usize,
    theta: f64,
) -> Result<Vec<[f64; 3]>, Error> {
    if missing == 0 {
        return Ok(Vec::new());
    }
    let axis = normalized(first)?;
    let reference = least_aligned_axis(axis);
    let first_perpendicular = normalized(cross(axis, reference))?;
    let second_perpendicular = cross(axis, first_perpendicular);
    let (sin_theta, cos_theta) = theta.sin_cos();
    let mut directions = Vec::with_capacity(missing);
    for step in 0..missing {
        let azimuth = TAU * step as f64 / missing as f64;
        let (sin_azimuth, cos_azimuth) = azimuth.sin_cos();
        let radial = add(
            scale(first_perpendicular, cos_azimuth),
            scale(second_perpendicular, sin_azimuth),
        );
        directions.push(normalized(add(
            scale(axis, cos_theta),
            scale(radial, sin_theta),
        ))?);
    }
    Ok(directions)
}

/// Directions completing the geometry around two existing bonds.
fn two_neighbour_directions(
    first: [f64; 3],
    second: [f64; 3],
    missing: usize,
    theta: f64,
) -> Result<Vec<[f64; 3]>, Error> {
    if missing == 0 {
        return Ok(Vec::new());
    }
    let first = normalized(first)?;
    let second = normalized(second)?;
    let bisector = normalized(add(first, second))?;
    let cos_alpha = dot(first, bisector).clamp(-1.0, 1.0);
    if cos_alpha.abs() <= DEGENERATE {
        return Err(Error::InvalidInput(
            "zero-length direction basis: the two bonds are perpendicular to their bisector"
                .to_owned(),
        ));
    }
    let alpha = cos_alpha.acos();
    let in_plane = normalized(add(first, scale(bisector, -cos_alpha)))?;
    match missing {
        1 => {
            let beta = PI - theta - alpha;
            if beta < 0.0 {
                return Ok(vec![scale(bisector, -1.0)]);
            }
            let (sin_beta, cos_beta) = beta.sin_cos();
            Ok(vec![normalized(add(
                scale(bisector, -cos_beta),
                scale(in_plane, sin_beta),
            ))?])
        }
        2 => {
            let cos_theta = theta.cos();
            let cos_beta = (-cos_theta / cos_alpha).clamp(-1.0, 1.0);
            let beta = cos_beta.acos();
            let normal = normalized(cross(bisector, in_plane))?;
            let (sin_beta, cos_beta) = beta.sin_cos();
            let along = scale(bisector, -cos_beta);
            Ok(vec![
                normalized(add(along, scale(normal, sin_beta)))?,
                normalized(add(along, scale(normal, -sin_beta)))?,
            ])
        }
        _ => Err(Error::InvalidInput(format!(
            "cannot place {missing} hydrogens on an atom with two neighbours"
        ))),
    }
}

/// Axis with the smallest component of `direction`, for a stable basis.
fn least_aligned_axis(direction: [f64; 3]) -> [f64; 3] {
    let [x, y, z] = direction.map(f64::abs);
    if x <= y && x <= z {
        [1.0, 0.0, 0.0]
    } else if y <= z {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    }
}

fn add(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[0] + second[0],
        first[1] + second[1],
        first[2] + second[2],
    ]
}

fn scale(vector: [f64; 3], factor: f64) -> [f64; 3] {
    [vector[0] * factor, vector[1] * factor, vector[2] * factor]
}

fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    first[0] * second[0] + first[1] * second[1] + first[2] * second[2]
}

fn cross(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    ]
}

pub(super) fn normalized(vector: [f64; 3]) -> Result<[f64; 3], Error> {
    let length = dot(vector, vector).sqrt();
    if !length.is_finite() || length <= DEGENERATE {
        return Err(Error::InvalidInput(
            "zero-length direction basis".to_owned(),
        ));
    }
    Ok(scale(vector, 1.0 / length))
}
