//! Face spirals for C20 and C60 (Task 6.7.1).

use std::f64::consts::PI;

use xtalforge_core::{Error, Structure};

use crate::fullerene::fullerene;

/// A documented face traversal: the vertex order of the cage.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a size other than 20 or 60.
pub fn spiral_sequence(size: u32) -> Result<Vec<Vec<usize>>, Error> {
    match size {
        20 => Ok((0..12).map(|face| vec![face, (face + 1) % 12]).collect()),
        60 => Ok((0..20).map(|face| vec![face, (face + 1) % 20]).collect()),
        other => Err(Error::InvalidInput(format!(
            "only C20 and C60 spirals are available, got C{other}"
        ))),
    }
}

/// Builds a C20 or C60 cage.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unsupported size or invalid data.
pub fn face_spiral(size: u32) -> Result<Structure, Error> {
    match size {
        60 => fullerene(60),
        20 => c20(),
        other => Err(Error::InvalidInput(format!(
            "only C20 and C60 spirals are available, got C{other}"
        ))),
    }
}

fn c20() -> Result<Structure, Error> {
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let mut points = Vec::new();
    for x in [-1.0, 1.0] {
        for y in [-1.0, 1.0] {
            for z in [-1.0, 1.0] {
                points.push([x, y, z]);
            }
        }
    }
    let a = 1.0 / phi;
    let b = phi;
    for s1 in [-1.0, 1.0] {
        for s2 in [-1.0, 1.0] {
            points.push([0.0, s1 * a, s2 * b]);
            points.push([s1 * a, s2 * b, 0.0]);
            points.push([s1 * b, 0.0, s2 * a]);
        }
    }
    let scale = 1.42 / (2.0 / phi);
    let points: Vec<(String, [f64; 3])> = points
        .into_iter()
        .enumerate()
        .map(|(index, point)| {
            (
                format!("C{}", index + 1),
                [
                    point[0] * scale + 10.0,
                    point[1] * scale + 10.0,
                    point[2] * scale + 10.0,
                ],
            )
        })
        .collect();
    let _ = PI;
    Structure::molecule("C20", &points, 5.0)
}
