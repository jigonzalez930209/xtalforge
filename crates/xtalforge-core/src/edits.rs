#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Structure edits: strain, rattle, site operations and optimal supercells.
//!
//! Every operation is additive and deterministic. The random operations take
//! an explicit seed (AGENTS.md rule 6) and use a small local SplitMix64
//! generator so the results do not depend on the platform.

use std::collections::BTreeMap;

use crate::{Error, Lattice, Site, Structure, StructureBuilder};

/// Converts engineering Voigt strain `[exx, eyy, ezz, gyz, gxz, gxy]` to a
/// symmetric tensor.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for non-finite values.
pub fn voigt_to_tensor(voigt: [f64; 6]) -> Result<[[f64; 3]; 3], Error> {
    if !voigt.iter().all(|value| value.is_finite()) {
        return Err(Error::InvalidInput(
            "strain components must be finite".to_owned(),
        ));
    }
    let [xx, yy, zz, yz, xz, xy] = voigt;
    Ok([
        [xx, 0.5 * xy, 0.5 * xz],
        [0.5 * xy, yy, 0.5 * yz],
        [0.5 * xz, 0.5 * yz, zz],
    ])
}

/// Converts a symmetric strain tensor to engineering Voigt components.
#[must_use]
pub fn tensor_to_voigt(tensor: [[f64; 3]; 3]) -> [f64; 6] {
    [
        tensor[0][0],
        tensor[1][1],
        tensor[2][2],
        2.0 * tensor[1][2],
        2.0 * tensor[0][2],
        2.0 * tensor[0][1],
    ]
}

/// Applies a Cartesian strain tensor to the cell; fractional coordinates stay.
///
/// # Errors
///
/// Returns the errors of [`Lattice::strained`] or a builder error.
pub fn strain(structure: &Structure, tensor: [[f64; 3]; 3]) -> Result<Structure, Error> {
    let lattice = structure.lattice().strained(tensor)?;
    copy_with_lattice(structure, lattice)
}

/// Applies a deformation gradient `F` (`x' = F x`).
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-finite or singular gradient, or a
/// builder error.
pub fn deform(structure: &Structure, gradient: [[f64; 3]; 3]) -> Result<Structure, Error> {
    if !gradient.iter().flatten().all(|value| value.is_finite()) {
        return Err(Error::InvalidInput(
            "deformation gradient must be finite".to_owned(),
        ));
    }
    let determinant = determinant(gradient);
    if determinant.abs() < 1e-9 || determinant <= 0.0 {
        return Err(Error::InvalidInput(
            "deformation gradient must have a positive determinant".to_owned(),
        ));
    }
    // `Lattice::strained` computes M (I + eps)^T; with eps = F^T - I the
    // result is M F, which is the row form of x' = F x.
    let transpose = [
        [gradient[0][0], gradient[1][0], gradient[2][0]],
        [gradient[0][1], gradient[1][1], gradient[2][1]],
        [gradient[0][2], gradient[1][2], gradient[2][2]],
    ];
    let eps = [
        [transpose[0][0] - 1.0, transpose[0][1], transpose[0][2]],
        [transpose[1][0], transpose[1][1] - 1.0, transpose[1][2]],
        [transpose[2][0], transpose[2][1], transpose[2][2] - 1.0],
    ];
    strain(structure, eps)
}

/// Isotropically scales the cell to a target volume in cubic Angstrom.
///
/// # Errors
///
/// Returns the errors of [`Lattice::scaled_to_volume`] or a builder error.
pub fn scale_to_volume(structure: &Structure, volume: f64) -> Result<Structure, Error> {
    let lattice = structure.lattice().scaled_to_volume(volume)?;
    copy_with_lattice(structure, lattice)
}

/// Isotropically scales the cell to a target density in g/cm3.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-positive density, or the errors
/// of [`Structure::cell_mass`] and [`scale_to_volume`].
pub fn scale_to_density(structure: &Structure, density: f64) -> Result<Structure, Error> {
    if !density.is_finite() || density <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "target density must be positive and finite, got {density}"
        )));
    }
    let mass = structure.cell_mass()?;
    let volume = mass / (density * 0.602_214_076_9);
    scale_to_volume(structure, volume)
}

/// Options for [`rattle`].
#[derive(Debug, Clone, PartialEq)]
pub struct RattleOptions {
    /// Gaussian standard deviation in Angstrom (or the maximum norm when
    /// `bounded` is set).
    pub amplitude: f64,
    /// Explicit seed for the deterministic generator.
    pub seed: u64,
    /// Draw a uniform direction with a fixed norm instead of a Gaussian.
    pub bounded: bool,
    /// Per-element amplitude overrides in Angstrom.
    pub per_element: BTreeMap<String, f64>,
}

impl Default for RattleOptions {
    fn default() -> Self {
        Self {
            amplitude: 0.05,
            seed: 0,
            bounded: false,
            per_element: BTreeMap::new(),
        }
    }
}

/// The rattled structure and its displacement statistics.
#[derive(Debug, Clone, PartialEq)]
pub struct RattleReport {
    /// The displaced P1 structure.
    pub structure: Structure,
    /// Largest Cartesian displacement in Angstrom.
    pub max_displacement: f64,
    /// Root mean square Cartesian displacement in Angstrom.
    pub rms_displacement: f64,
    /// Number of displaced atoms.
    pub n_atoms: usize,
}

/// Randomly displaces every atom of the P1 cell.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-finite amplitude, or a builder
/// error.
pub fn rattle(structure: &Structure, options: &RattleOptions) -> Result<RattleReport, Error> {
    if !options.amplitude.is_finite() || options.amplitude < 0.0 {
        return Err(Error::InvalidInput(
            "rattle amplitude must be non-negative and finite".to_owned(),
        ));
    }
    for (element, amplitude) in &options.per_element {
        if !amplitude.is_finite() || *amplitude < 0.0 {
            return Err(Error::InvalidInput(format!(
                "rattle amplitude for {element} must be non-negative and finite"
            )));
        }
    }
    let p1 = structure.p1()?;
    let mut rng = SplitMix64::new(options.seed);
    let mut builder = StructureBuilder::new(structure.name(), p1.lattice);
    let mut squared = 0.0;
    let mut largest = 0.0_f64;
    for atom in &p1.atoms {
        let amplitude = options
            .per_element
            .get(&atom.element)
            .copied()
            .unwrap_or(options.amplitude);
        let displacement = if options.bounded {
            let direction = rng.unit_vector();
            [
                direction[0] * amplitude,
                direction[1] * amplitude,
                direction[2] * amplitude,
            ]
        } else {
            let scale = amplitude * (rng.gaussian() / std::f64::consts::SQRT_2);
            let direction = rng.unit_vector();
            [
                direction[0] * scale,
                direction[1] * scale,
                direction[2] * scale,
            ]
        };
        let cart = p1.lattice.to_cart(atom.frac);
        let moved = [
            cart[0] + displacement[0],
            cart[1] + displacement[1],
            cart[2] + displacement[2],
        ];
        let frac = Lattice::wrap(p1.lattice.to_frac(moved));
        let site =
            Site::new(atom.label.clone(), &atom.element, frac)?.with_occupancy(atom.occupancy)?;
        builder = builder.site(site);
        let norm = (displacement[0] * displacement[0]
            + displacement[1] * displacement[1]
            + displacement[2] * displacement[2])
            .sqrt();
        squared += norm * norm;
        largest = largest.max(norm);
    }
    Ok(RattleReport {
        structure: builder.build()?,
        max_displacement: largest,
        rms_displacement: if p1.atoms.is_empty() {
            0.0
        } else {
            (squared / p1.atoms.len() as f64).sqrt()
        },
        n_atoms: p1.atoms.len(),
    })
}

mod sites;
mod supercell;

pub use sites::{
    SiteSort, SpeciesReplacement, add_sites, move_site, remove_sites, replace_site,
    replace_species, sort_sites, wrap_sites,
};
pub use supercell::{SupercellChoice, optimal_supercell};

fn copy_with_lattice(structure: &Structure, lattice: Lattice) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new(structure.name(), lattice).hall(structure.hall());
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for site in structure.sites() {
        builder = builder.site(site.clone());
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

/// Deterministic SplitMix64 generator.
#[derive(Debug, Clone)]
pub(crate) struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    pub(crate) fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    pub(crate) fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    pub(crate) fn gaussian(&mut self) -> f64 {
        let u1 = self.uniform().max(1e-300);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }

    pub(crate) fn unit_vector(&mut self) -> [f64; 3] {
        let z = 2.0 * self.uniform() - 1.0;
        let phi = 2.0 * std::f64::consts::PI * self.uniform();
        let radius = (1.0 - z * z).max(0.0).sqrt();
        [radius * phi.cos(), radius * phi.sin(), z]
    }

    pub(crate) fn shuffle<T>(&mut self, values: &mut [T]) {
        for index in (1..values.len()).rev() {
            let other = (self.next_u64() % (index as u64 + 1)) as usize;
            values.swap(index, other);
        }
    }
}

fn determinant(matrix: [[f64; 3]; 3]) -> f64 {
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
}
