//! Distances, angles, torsions, centroids and planes (Task 1.14.1).

use crate::{Error, Structure};

/// Minimum-image distance between two atoms in the P1 cell.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an out-of-range index or a failed P1
/// expansion.
pub fn distance(structure: &Structure, i: usize, j: usize) -> Result<f64, Error> {
    let p1 = structure.p1()?;
    let first = p1
        .atoms
        .get(i)
        .ok_or_else(|| index_error(i, p1.atoms.len()))?;
    let second = p1
        .atoms
        .get(j)
        .ok_or_else(|| index_error(j, p1.atoms.len()))?;
    Ok(p1.lattice.distance(first.frac, second.frac))
}

/// Angle `i-j-k` in degrees, using minimum images around `j`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an out-of-range index or a degenerate
/// angle.
pub fn angle(structure: &Structure, i: usize, j: usize, k: usize) -> Result<f64, Error> {
    let p1 = structure.p1()?;
    let pivot = p1
        .atoms
        .get(j)
        .ok_or_else(|| index_error(j, p1.atoms.len()))?;
    let first = p1
        .atoms
        .get(i)
        .ok_or_else(|| index_error(i, p1.atoms.len()))?;
    let second = p1
        .atoms
        .get(k)
        .ok_or_else(|| index_error(k, p1.atoms.len()))?;
    let first_vector = p1
        .lattice
        .to_cart(p1.lattice.minimum_image(subtract(first.frac, pivot.frac)).0);
    let second_vector = p1.lattice.to_cart(
        p1.lattice
            .minimum_image(subtract(second.frac, pivot.frac))
            .0,
    );
    let first_norm = norm(first_vector);
    let second_norm = norm(second_vector);
    if first_norm < 1e-12 || second_norm < 1e-12 {
        return Err(Error::InvalidInput(
            "angle has a zero-length arm".to_owned(),
        ));
    }
    let cosine = dot(first_vector, second_vector) / (first_norm * second_norm);
    Ok(cosine.clamp(-1.0, 1.0).acos().to_degrees())
}

/// Torsion `i-j-k-l` in degrees.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an out-of-range index or a degenerate
/// torsion.
pub fn torsion(
    structure: &Structure,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
) -> Result<f64, Error> {
    let p1 = structure.p1()?;
    let point = |index: usize| {
        p1.atoms
            .get(index)
            .map(|atom| atom.frac)
            .ok_or_else(|| index_error(index, p1.atoms.len()))
    };
    let fi = point(i)?;
    let fj = point(j)?;
    let fk = point(k)?;
    let fl = point(l)?;
    let b1 = p1
        .lattice
        .to_cart(p1.lattice.minimum_image(subtract(fj, fi)).0);
    let b2 = p1
        .lattice
        .to_cart(p1.lattice.minimum_image(subtract(fk, fj)).0);
    let b3 = p1
        .lattice
        .to_cart(p1.lattice.minimum_image(subtract(fl, fk)).0);
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let b2_norm = norm(b2);
    if norm(n1) < 1e-12 || norm(n2) < 1e-12 || b2_norm < 1e-12 {
        return Err(Error::InvalidInput(
            "torsion has a collinear segment".to_owned(),
        ));
    }
    let m = cross(n1, scale(b2, 1.0 / b2_norm));
    let x = dot(n1, n2);
    let y = dot(m, n2);
    Ok(y.atan2(x).to_degrees())
}

/// Centroid of Cartesian points.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for no points.
pub fn centroid(points: &[[f64; 3]]) -> Result<[f64; 3], Error> {
    if points.is_empty() {
        return Err(Error::InvalidInput(
            "a centroid needs at least one point".to_owned(),
        ));
    }
    let mut sum = [0.0; 3];
    for point in points {
        sum[0] += point[0];
        sum[1] += point[1];
        sum[2] += point[2];
    }
    let count = points.len() as f64;
    Ok([sum[0] / count, sum[1] / count, sum[2] / count])
}

/// Angle between two planes in degrees.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for zero normals.
pub fn interplanar_angle(first_normal: [f64; 3], second_normal: [f64; 3]) -> Result<f64, Error> {
    let first_norm = norm(first_normal);
    let second_norm = norm(second_normal);
    if first_norm < 1e-12 || second_norm < 1e-12 {
        return Err(Error::InvalidInput(
            "a plane normal has zero length".to_owned(),
        ));
    }
    let cosine = dot(first_normal, second_normal) / (first_norm * second_norm);
    Ok(cosine.clamp(-1.0, 1.0).acos().to_degrees())
}

fn index_error(index: usize, length: usize) -> Error {
    Error::InvalidInput(format!(
        "atom index {index} is outside the P1 cell (length {length})"
    ))
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
