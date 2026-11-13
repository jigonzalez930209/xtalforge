//! Bravais classification and band paths (Task 2.10.2-2.10.4).
//!
//! The special points and path segments are generated from ASE 3.29 and stored
//! in `data/kpaths.json.gz`; the table is decompressed and parsed on first use.
//! Coordinates are fractional in the reciprocal basis of the user's cell and
//! are converted to Cartesian `1/Angstrom` (with `2 pi`) on request.

use std::collections::BTreeMap;

use serde::Deserialize;
use xtalforge_core::data::{EmbeddedBundle, JsonBundle};
use xtalforge_core::{Error, Lattice, Structure};

use crate::dataset::{SymmetryOptions, detect_symmetry};
use crate::spacegroup::SpaceGroup;

/// One Bravais lattice entry.
#[derive(Debug, Clone, Deserialize)]
struct Entry {
    special_points: BTreeMap<String, [f64; 3]>,
    path: Vec<Vec<String>>,
}

#[derive(Debug, Deserialize)]
struct Database {
    lattices: BTreeMap<String, Entry>,
}

/// The gzip path table is decompressed, checksum-verified and parsed on first
/// use (Task 15.3.2); earlier code parsed it on every call.
static DATABASE: JsonBundle<Database> = JsonBundle::new(EmbeddedBundle::new(
    "xtalforge-symmetry/kpaths.json.gz",
    include_bytes!("../data/kpaths.json.gz"),
    "c36354ee1946693dc8bcdd587b296b9a96418ca84a6dad64d60b04f596b85676",
));

/// A band path in the user's reciprocal basis.
#[derive(Debug, Clone, PartialEq)]
pub struct KPath {
    /// Bravais label, for example `cF`.
    pub bravais: String,
    /// Special points as fractional reciprocal coordinates.
    pub special_points: BTreeMap<String, [f64; 3]>,
    /// Path segments as label pairs.
    pub path: Vec<(String, String)>,
}

/// Classifies a lattice into one of the 14 Bravais types.
///
/// The classification uses the conventional metric and the primitive volume
/// ratio from the Niggli reduction to distinguish centering.
#[must_use]
pub fn classify_bravais(lattice: &Lattice) -> &'static str {
    let [a, b, c, alpha, beta, gamma] = lattice.parameters();
    let angle_90 = |value: f64| (value - 90.0).abs() < 1e-3;
    let equal = |first: f64, second: f64| (first - second).abs() < 1e-3;
    let ratio = primitive_ratio(lattice);
    if equal(a, b) && equal(b, c) {
        if angle_90(alpha) && angle_90(beta) && angle_90(gamma) {
            return if ratio > 3.0 {
                "cF"
            } else if ratio > 1.5 {
                "cI"
            } else {
                "cP"
            };
        }
        if equal(alpha, beta) && equal(beta, gamma) {
            if (alpha - 60.0).abs() < 1e-3 {
                return "cF";
            }
            if (alpha - 109.471_220_634_490_69).abs() < 1e-3 {
                return "cI";
            }
            return "hR";
        }
    }
    if equal(a, b) && angle_90(alpha) && angle_90(beta) && (gamma - 120.0).abs() < 1e-3 {
        return "hP";
    }
    if angle_90(alpha) && angle_90(beta) && angle_90(gamma) {
        if equal(a, b) && !equal(b, c) {
            return if ratio > 1.5 { "tI" } else { "tP" };
        }
        if !equal(a, b) && !equal(b, c) {
            return if ratio > 3.0 {
                "oF"
            } else if ratio > 1.5 {
                "oC"
            } else {
                "oP"
            };
        }
    }
    let right = [alpha, beta, gamma]
        .iter()
        .filter(|angle| angle_90(**angle))
        .count();
    if right == 2 { "mP" } else { "aP" }
}

/// Returns the band path for a lattice alone.
///
/// The centering cannot be inferred from a conventional lattice without the
/// basis; use [`band_path_for_structure`] when a structure is available.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the Bravais entry is missing.
pub fn band_path(lattice: &Lattice) -> Result<KPath, Error> {
    band_path_for_bravais(classify_bravais(lattice))
}

/// Returns the band path for a structure, using the detected space group to
/// determine the Bravais centering.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the entry is missing or the errors of
/// [`detect_symmetry`].
pub fn band_path_for_structure(
    structure: &Structure,
    options: SymmetryOptions,
) -> Result<KPath, Error> {
    let dataset = detect_symmetry(structure, options)?;
    let group = SpaceGroup::from_hall_number(dataset.hall_number)?;
    band_path_for_bravais(&group.bravais)
}

fn band_path_for_bravais(bravais: &str) -> Result<KPath, Error> {
    let database = DATABASE.get()?;
    let entry = database.lattices.get(bravais).ok_or_else(|| {
        Error::InvalidInput(format!("no k-path entry for Bravais type {bravais}"))
    })?;
    Ok(KPath {
        bravais: bravais.to_owned(),
        special_points: entry.special_points.clone(),
        path: entry
            .path
            .iter()
            .flat_map(|segment| {
                segment
                    .windows(2)
                    .filter_map(|pair| Some((pair.first()?.clone(), pair.get(1)?.clone())))
            })
            .collect(),
    })
}

impl KPath {
    /// Cartesian `k` vectors of the special points in `1/Angstrom` with `2 pi`.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`Lattice::reciprocal`].
    pub fn cartesian_points(&self, lattice: &Lattice) -> Result<BTreeMap<String, [f64; 3]>, Error> {
        let reciprocal = lattice.reciprocal(true)?;
        Ok(self
            .special_points
            .keys()
            .cloned()
            .zip(
                self.special_points
                    .values()
                    .map(|point| reciprocal.to_cart(*point)),
            )
            .collect())
    }

    /// Vertices along the path with the cumulative Cartesian distance.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`KPath::cartesian_points`].
    pub fn vertices(&self, lattice: &Lattice) -> Result<Vec<(String, [f64; 3], f64)>, Error> {
        let points = self.cartesian_points(lattice)?;
        let mut vertices = Vec::new();
        let mut distance = 0.0;
        let mut previous: Option<[f64; 3]> = None;
        for (start, end) in &self.path {
            for label in [start, end] {
                let Some(point) = points.get(label).copied() else {
                    continue;
                };
                if let Some(previous) = previous {
                    distance += norm(subtract(point, previous));
                }
                vertices.push((label.clone(), point, distance));
                previous = Some(point);
            }
        }
        Ok(vertices)
    }
}

fn primitive_ratio(lattice: &Lattice) -> f64 {
    let moyo_lattice = moyo::base::Lattice::from_basis(lattice.matrix());
    let (reduced, _) = moyo_lattice.unchecked_niggli_reduce();
    let primitive_volume = reduced.basis.determinant().abs();
    if primitive_volume <= 0.0 {
        return 1.0;
    }
    lattice.volume() / primitive_volume
}

fn norm(vector: [f64; 3]) -> f64 {
    (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]).sqrt()
}

fn subtract(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[0] - second[0],
        first[1] - second[1],
        first[2] - second[2],
    ]
}
