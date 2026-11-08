//! Internal coordinates for scans and constraints (Task 4.5.1).
//!
//! Three coordinates are supported at this core level: a distance, an angle
//! and a torsion, each addressed by asymmetric-site indices.  Scans operate on
//! P1 structures, so those indices are also P1 atom indices.  Periodic images
//! are taken into account through the lattice's minimum-image convention.
//!
//! Values are Angstrom for a distance and degrees for an angle or torsion.
//! The analytic gradients are differentiated with respect to Cartesian site
//! positions; a degenerate angle or torsion has a zero gradient and is refused
//! explicitly when it is used as a scan axis.

use serde::{Deserialize, Serialize};
use xtalforge_core::Lattice;

use crate::Error;
use crate::math::{add, cross, dot, norm, scale, sub, unit};

/// One internal coordinate.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum ScanCoordinate {
    /// Distance between sites `i` and `j`.
    Distance {
        /// First site index.
        i: usize,
        /// Second site index.
        j: usize,
    },
    /// Angle at site `j` between `i` and `k`, in degrees.
    Angle {
        /// First endpoint.
        i: usize,
        /// Vertex.
        j: usize,
        /// Second endpoint.
        k: usize,
    },
    /// Torsion across `i-j-k-l`, in degrees in `(-180, 180]`.
    Torsion {
        /// First atom.
        i: usize,
        /// Second atom.
        j: usize,
        /// Third atom.
        k: usize,
        /// Fourth atom.
        l: usize,
    },
}

impl ScanCoordinate {
    /// Distance coordinate.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the two indices are equal.
    pub fn distance(i: usize, j: usize) -> Result<Self, Error> {
        if i == j {
            return Err(Error::InvalidInput(
                "a distance needs two different atoms".to_owned(),
            ));
        }
        Ok(Self::Distance { i, j })
    }

    /// Angle coordinate, vertex in the middle.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when two indices are equal.
    pub fn angle(i: usize, j: usize, k: usize) -> Result<Self, Error> {
        if i == j || j == k || i == k {
            return Err(Error::InvalidInput(
                "an angle needs three different atoms".to_owned(),
            ));
        }
        Ok(Self::Angle { i, j, k })
    }

    /// Torsion coordinate.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when two indices are equal.
    pub fn torsion(i: usize, j: usize, k: usize, l: usize) -> Result<Self, Error> {
        let indices = [i, j, k, l];
        for (first, value) in indices.iter().enumerate() {
            for other in indices.iter().skip(first + 1) {
                if value == other {
                    return Err(Error::InvalidInput(
                        "a torsion needs four different atoms".to_owned(),
                    ));
                }
            }
        }
        Ok(Self::Torsion { i, j, k, l })
    }

    /// Human-readable label.
    #[must_use]
    pub fn label(&self) -> String {
        match self {
            Self::Distance { i, j } => format!("distance {i}-{j}"),
            Self::Angle { i, j, k } => format!("angle {i}-{j}-{k}"),
            Self::Torsion { i, j, k, l } => format!("torsion {i}-{j}-{k}-{l}"),
        }
    }

    /// Units of the value: `"A"` or `"deg"`.
    #[must_use]
    pub fn units(&self) -> &'static str {
        match self {
            Self::Distance { .. } => "A",
            Self::Angle { .. } | Self::Torsion { .. } => "deg",
        }
    }

    /// Site indices this coordinate reads.
    #[must_use]
    pub fn indices(&self) -> Vec<usize> {
        match self {
            Self::Distance { i, j } => vec![*i, *j],
            Self::Angle { i, j, k } => vec![*i, *j, *k],
            Self::Torsion { i, j, k, l } => vec![*i, *j, *k, *l],
        }
    }

    /// Checks the indices against a site count.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an index outside `0..n_sites`.
    pub fn validate(&self, n_sites: usize) -> Result<(), Error> {
        for index in self.indices() {
            if index >= n_sites {
                return Err(Error::InvalidInput(format!(
                    "{} refers to atom {index}, but the structure has {n_sites} sites",
                    self.label()
                )));
            }
        }
        Ok(())
    }

    /// Coordinate value for Cartesian site positions.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an index outside the position
    /// array, a degenerate geometry, or a non-finite result.
    pub fn value(&self, positions: &[[f64; 3]], lattice: &Lattice) -> Result<f64, Error> {
        match self {
            Self::Distance { i, j } => {
                let vector = minimum_vector(positions, lattice, *i, *j)?;
                let distance = norm(vector);
                if distance <= 1e-12 {
                    return Err(Error::InvalidInput(
                        "the two atoms of a distance coordinate coincide".to_owned(),
                    ));
                }
                Ok(distance)
            }
            Self::Angle { i, j, k } => {
                let u = minimum_vector(positions, lattice, *j, *i)?;
                let v = minimum_vector(positions, lattice, *j, *k)?;
                let denominator = norm(u) * norm(v);
                if denominator <= 1e-14 {
                    return Err(Error::InvalidInput(
                        "an angle coordinate has a zero-length arm".to_owned(),
                    ));
                }
                let cosine = (dot(u, v) / denominator).clamp(-1.0, 1.0);
                Ok(cosine.acos().to_degrees())
            }
            Self::Torsion { i, j, k, l } => {
                let b1 = minimum_vector(positions, lattice, *i, *j)?;
                let b2 = minimum_vector(positions, lattice, *j, *k)?;
                let b3 = minimum_vector(positions, lattice, *k, *l)?;
                let n1 = cross(b1, b2);
                let n2 = cross(b2, b3);
                let length = norm(b2);
                if length < 1e-9 || norm(n1) < 1e-9 || norm(n2) < 1e-9 {
                    return Err(Error::InvalidInput(
                        "a torsion coordinate is undefined where three atoms are collinear"
                            .to_owned(),
                    ));
                }
                let x = dot(n1, n2);
                let y = dot(cross(n1, n2), b2) / length;
                Ok(y.atan2(x).to_degrees())
            }
        }
    }

    /// Cartesian gradient `d(value)/d(position)` for every site.
    ///
    /// # Errors
    ///
    /// Returns the same errors as [`ScanCoordinate::value`].
    pub fn gradient(
        &self,
        positions: &[[f64; 3]],
        lattice: &Lattice,
    ) -> Result<Vec<[f64; 3]>, Error> {
        let mut result = vec![[0.0; 3]; positions.len()];
        match self {
            Self::Distance { i, j } => {
                let vector = minimum_vector(positions, lattice, *i, *j)?;
                let Some(direction) = unit(vector) else {
                    return Err(Error::InvalidInput(
                        "the two atoms of a distance coordinate coincide".to_owned(),
                    ));
                };
                scatter(&mut result, *i, scale(direction, -1.0));
                scatter(&mut result, *j, direction);
            }
            Self::Angle { i, j, k } => {
                let u = minimum_vector(positions, lattice, *j, *i)?;
                let v = minimum_vector(positions, lattice, *j, *k)?;
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
                scatter(&mut result, *i, first);
                scatter(&mut result, *j, scale(add(first, third), -1.0));
                scatter(&mut result, *k, third);
            }
            Self::Torsion { i, j, k, l } => {
                let b1 = minimum_vector(positions, lattice, *i, *j)?;
                let b2 = minimum_vector(positions, lattice, *j, *k)?;
                let b3 = minimum_vector(positions, lattice, *k, *l)?;
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
                scatter(&mut result, *i, first);
                scatter(&mut result, *j, second);
                scatter(&mut result, *k, third);
                scatter(&mut result, *l, fourth);
            }
        }
        Ok(result)
    }

    /// Moves the geometry so the coordinate reads `target`.
    ///
    /// This is a seed, not a constraint: it translates or rotates one atom to
    /// the requested value and leaves the rest for the optimizer.  For a
    /// torsion it rotates the last atom about the central bond; for an angle
    /// it rotates the third atom about the normal at the vertex.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a degenerate starting geometry or
    /// a non-finite target.
    pub fn set(
        &self,
        positions: &mut [[f64; 3]],
        lattice: &Lattice,
        target: f64,
    ) -> Result<(), Error> {
        if !target.is_finite() {
            return Err(Error::InvalidInput("scan target must be finite".to_owned()));
        }
        match self {
            Self::Distance { i, j } => {
                let vector = minimum_vector(positions, lattice, *i, *j)?;
                let Some(direction) = unit(vector) else {
                    return Err(Error::InvalidInput(
                        "the two atoms of a distance coordinate coincide".to_owned(),
                    ));
                };
                let change = target - norm(vector);
                if let Some(position) = positions.get_mut(*j) {
                    *position = add(*position, scale(direction, change));
                }
            }
            Self::Angle { i, j, k } => {
                let u = minimum_vector(positions, lattice, *j, *i)?;
                let v = minimum_vector(positions, lattice, *j, *k)?;
                let Some(axis) = unit(cross(u, v)) else {
                    return Err(Error::InvalidInput(
                        "cannot set an angle whose three atoms are collinear".to_owned(),
                    ));
                };
                let current = self.value(positions, lattice)?;
                let rotated = rotate(v, axis, (target - current).to_radians());
                let vertex = positions.get(*j).copied().unwrap_or([0.0; 3]);
                if let Some(position) = positions.get_mut(*k) {
                    *position = add(vertex, rotated);
                }
            }
            Self::Torsion { j, k, l, .. } => {
                let b2 = minimum_vector(positions, lattice, *j, *k)?;
                let Some(axis) = unit(b2) else {
                    return Err(Error::InvalidInput(
                        "cannot set a torsion whose central bond is zero".to_owned(),
                    ));
                };
                let pivot = positions.get(*k).copied().unwrap_or([0.0; 3]);
                let arm = minimum_vector(positions, lattice, *k, *l)?;
                let current = self.value(positions, lattice)?;
                let rotated = rotate(arm, axis, (target - current).to_radians());
                if let Some(position) = positions.get_mut(*l) {
                    *position = add(pivot, rotated);
                }
            }
        }
        Ok(())
    }
}

/// Minimum-image Cartesian vector from site `first` to site `second`.
fn minimum_vector(
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

/// Adds a contribution to one site's gradient.
fn scatter(gradient: &mut [[f64; 3]], index: usize, value: [f64; 3]) {
    if let Some(slot) = gradient.get_mut(index) {
        *slot = add(*slot, value);
    }
}

/// Rodrigues rotation of `vector` about a unit `axis`.
fn rotate(vector: [f64; 3], axis: [f64; 3], angle: f64) -> [f64; 3] {
    let cosine = angle.cos();
    let sine = angle.sin();
    let parallel = scale(axis, dot(axis, vector) * (1.0 - cosine));
    let perpendicular = scale(cross(axis, vector), sine);
    add(add(scale(vector, cosine), parallel), perpendicular)
}
