//! Distance, angle and torsion value, gradient and seeding.
//!
//! These are the original Task 4.5 coordinates.  They live here so
//! [`super::coordinates`] can stay a dispatch layer that also covers cell
//! parameters and collective variables (Task 14.2.3).
//!
//! Values are Angstrom for a distance and degrees for an angle or torsion.
//! Gradients are with respect to Cartesian site positions, and a degenerate
//! geometry returns a zero gradient and is refused when used as a scan axis.

use xtalforge_core::Lattice;

use crate::Error;
use crate::math::{add, cross, dot, norm, scale, sub, unit};

/// Distance between sites `i` and `j`, in Angstrom.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an out-of-range index or coincident
/// atoms.
pub(crate) fn distance_value(
    positions: &[[f64; 3]],
    lattice: &Lattice,
    i: usize,
    j: usize,
) -> Result<f64, Error> {
    let vector = minimum_vector(positions, lattice, i, j)?;
    let distance = norm(vector);
    if distance <= 1e-12 {
        return Err(Error::InvalidInput(
            "the two atoms of a distance coordinate coincide".to_owned(),
        ));
    }
    Ok(distance)
}

/// Angle at `j` between `i` and `k`, in degrees.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an out-of-range index or a zero arm.
pub(crate) fn angle_value(
    positions: &[[f64; 3]],
    lattice: &Lattice,
    i: usize,
    j: usize,
    k: usize,
) -> Result<f64, Error> {
    let u = minimum_vector(positions, lattice, j, i)?;
    let v = minimum_vector(positions, lattice, j, k)?;
    let denominator = norm(u) * norm(v);
    if denominator <= 1e-14 {
        return Err(Error::InvalidInput(
            "an angle coordinate has a zero-length arm".to_owned(),
        ));
    }
    let cosine = (dot(u, v) / denominator).clamp(-1.0, 1.0);
    Ok(cosine.acos().to_degrees())
}

/// Torsion across `i-j-k-l`, in degrees in `(-180, 180]`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an out-of-range index or collinear
/// atoms.
pub(crate) fn torsion_value(
    positions: &[[f64; 3]],
    lattice: &Lattice,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
) -> Result<f64, Error> {
    let (b2, n1, n2) = torsion_frames(positions, lattice, i, j, k, l)?;
    let length = norm(b2);
    if length < 1e-9 || norm(n1) < 1e-9 || norm(n2) < 1e-9 {
        return Err(Error::InvalidInput(
            "a torsion coordinate is undefined where three atoms are collinear".to_owned(),
        ));
    }
    let x = dot(n1, n2);
    let y = dot(cross(n1, n2), b2) / length;
    Ok(y.atan2(x).to_degrees())
}

/// Cartesian gradient of the distance coordinate.
pub(crate) fn distance_gradient(
    positions: &[[f64; 3]],
    lattice: &Lattice,
    i: usize,
    j: usize,
) -> Result<Vec<[f64; 3]>, Error> {
    let mut result = vec![[0.0; 3]; positions.len()];
    let vector = minimum_vector(positions, lattice, i, j)?;
    let Some(direction) = unit(vector) else {
        return Err(Error::InvalidInput(
            "the two atoms of a distance coordinate coincide".to_owned(),
        ));
    };
    scatter(&mut result, i, scale(direction, -1.0));
    scatter(&mut result, j, direction);
    Ok(result)
}

/// Cartesian gradient of the angle coordinate.
pub(crate) fn angle_gradient(
    positions: &[[f64; 3]],
    lattice: &Lattice,
    i: usize,
    j: usize,
    k: usize,
) -> Result<Vec<[f64; 3]>, Error> {
    let mut result = vec![[0.0; 3]; positions.len()];
    let u = minimum_vector(positions, lattice, j, i)?;
    let v = minimum_vector(positions, lattice, j, k)?;
    let lu = norm(u);
    let lv = norm(v);
    if lu < 1e-12 || lv < 1e-12 {
        return Err(Error::InvalidInput(
            "an angle coordinate has a zero-length arm".to_owned(),
        ));
    }
    let uh = scale(u, 1.0 / lu);
    let vh = scale(v, 1.0 / lv);
    let cosine = dot(uh, vh).clamp(-1.0, 1.0);
    let sine = (1.0 - cosine * cosine).sqrt();
    if sine < 1e-9 {
        return Ok(result);
    }
    let factor = 180.0 / std::f64::consts::PI;
    let first = scale(sub(scale(uh, cosine), vh), factor / (lu * sine));
    let third = scale(sub(scale(vh, cosine), uh), factor / (lv * sine));
    scatter(&mut result, i, first);
    scatter(&mut result, j, scale(add(first, third), -1.0));
    scatter(&mut result, k, third);
    Ok(result)
}

/// Cartesian gradient of the torsion coordinate.
pub(crate) fn torsion_gradient(
    positions: &[[f64; 3]],
    lattice: &Lattice,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
) -> Result<Vec<[f64; 3]>, Error> {
    let mut result = vec![[0.0; 3]; positions.len()];
    let (b1, b2, b3) = torsion_bonds(positions, lattice, i, j, k, l)?;
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let l2 = norm(b2);
    let s1 = dot(n1, n1);
    let s2 = dot(n2, n2);
    if l2 < 1e-9 || s1 < 1e-18 || s2 < 1e-18 {
        return Ok(result);
    }
    let factor = 180.0 / std::f64::consts::PI;
    let first = scale(n1, -l2 / s1 * factor);
    let fourth = scale(n2, l2 / s2 * factor);
    let a = dot(b1, b2) / (l2 * l2);
    let b = dot(b3, b2) / (l2 * l2);
    let second = add(scale(first, -(a + 1.0)), scale(fourth, b));
    let third = sub(scale(first, a), scale(fourth, b + 1.0));
    scatter(&mut result, i, first);
    scatter(&mut result, j, second);
    scatter(&mut result, k, third);
    scatter(&mut result, l, fourth);
    Ok(result)
}

/// Moves the geometry so the distance reads `target`.
pub(crate) fn set_distance(
    positions: &mut [[f64; 3]],
    lattice: &Lattice,
    i: usize,
    j: usize,
    target: f64,
) -> Result<(), Error> {
    let vector = minimum_vector(positions, lattice, i, j)?;
    let Some(direction) = unit(vector) else {
        return Err(Error::InvalidInput(
            "the two atoms of a distance coordinate coincide".to_owned(),
        ));
    };
    let change = target - norm(vector);
    if let Some(position) = positions.get_mut(j) {
        *position = add(*position, scale(direction, change));
    }
    Ok(())
}

/// Moves the geometry so the angle reads `target`.
pub(crate) fn set_angle(
    positions: &mut [[f64; 3]],
    lattice: &Lattice,
    i: usize,
    j: usize,
    k: usize,
    target: f64,
) -> Result<(), Error> {
    let u = minimum_vector(positions, lattice, j, i)?;
    let v = minimum_vector(positions, lattice, j, k)?;
    let Some(axis) = unit(cross(u, v)) else {
        return Err(Error::InvalidInput(
            "cannot set an angle whose three atoms are collinear".to_owned(),
        ));
    };
    let current = angle_value(positions, lattice, i, j, k)?;
    let rotated = rotate(v, axis, (target - current).to_radians());
    let vertex = positions.get(j).copied().unwrap_or([0.0; 3]);
    if let Some(position) = positions.get_mut(k) {
        *position = add(vertex, rotated);
    }
    Ok(())
}

/// Moves the geometry so the torsion reads `target`.
pub(crate) fn set_torsion(
    positions: &mut [[f64; 3]],
    lattice: &Lattice,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
    target: f64,
) -> Result<(), Error> {
    let b2 = minimum_vector(positions, lattice, j, k)?;
    let Some(axis) = unit(b2) else {
        return Err(Error::InvalidInput(
            "cannot set a torsion whose central bond is zero".to_owned(),
        ));
    };
    let pivot = positions.get(k).copied().unwrap_or([0.0; 3]);
    let arm = minimum_vector(positions, lattice, k, l)?;
    let current = torsion_value(positions, lattice, i, j, k, l)?;
    let rotated = rotate(arm, axis, (target - current).to_radians());
    if let Some(position) = positions.get_mut(l) {
        *position = add(pivot, rotated);
    }
    Ok(())
}

/// Minimum-image Cartesian vector from site `first` to site `second`.
pub(crate) fn minimum_vector(
    positions: &[[f64; 3]],
    lattice: &Lattice,
    first: usize,
    second: usize,
) -> Result<[f64; 3], Error> {
    let Some(from) = positions.get(first) else {
        return Err(Error::InvalidInput(format!(
            "atom index {first} is outside the position list"
        )));
    };
    let Some(to) = positions.get(second) else {
        return Err(Error::InvalidInput(format!(
            "atom index {second} is outside the position list"
        )));
    };
    let delta = sub(*to, *from);
    let frac = lattice.to_frac(delta);
    let (wrapped, _) = lattice.minimum_image(frac);
    Ok(lattice.to_cart(wrapped))
}

/// The three bonded vectors of a torsion: `i-j`, `j-k` and `k-l`.
type TorsionBonds = ([f64; 3], [f64; 3], [f64; 3]);

/// Central bond and the two plane normals of a torsion.
fn torsion_frames(
    positions: &[[f64; 3]],
    lattice: &Lattice,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
) -> Result<TorsionBonds, Error> {
    let (b1, b2, b3) = torsion_bonds(positions, lattice, i, j, k, l)?;
    Ok((b2, cross(b1, b2), cross(b2, b3)))
}

/// The three bonded vectors of a torsion.
fn torsion_bonds(
    positions: &[[f64; 3]],
    lattice: &Lattice,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
) -> Result<TorsionBonds, Error> {
    Ok((
        minimum_vector(positions, lattice, i, j)?,
        minimum_vector(positions, lattice, j, k)?,
        minimum_vector(positions, lattice, k, l)?,
    ))
}

/// Adds a contribution to one site's gradient.
fn scatter(gradient: &mut [[f64; 3]], index: usize, value: [f64; 3]) {
    if let Some(slot) = gradient.get_mut(index) {
        *slot = add(*slot, value);
    }
}

/// Rodrigues rotation of `vector` about a unit `axis`.
pub(crate) fn rotate(vector: [f64; 3], axis: [f64; 3], angle: f64) -> [f64; 3] {
    let cosine = angle.cos();
    let sine = angle.sin();
    let parallel = scale(axis, dot(axis, vector) * (1.0 - cosine));
    let perpendicular = scale(cross(axis, vector), sine);
    add(add(scale(vector, cosine), parallel), perpendicular)
}
