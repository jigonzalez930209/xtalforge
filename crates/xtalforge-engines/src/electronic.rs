//! Electronic-structure analysis (Task 18.3).
//!
//! The module consumes plain arrays through [`ElectronicStructure`] and returns
//! typed, serializable results. It does not depend on the parser result types
//! from Task 18.2: a parser, a recorded fixture or a test can build the input
//! directly.
//!
//! - [`band_structure`] turns eigenvalues on a k-path into bands and
//!   [`band_gap`] reports the direct/indirect gap; [`spin_gaps`] resolves both
//!   spin channels.
//! - [`effective_mass`] fits a parabola around a band edge and returns the mass
//!   in electron masses.
//! - [`density_of_states`] builds total, element- and orbital-projected DOS on
//!   an energy grid; [`smear_dos`] and [`d_band_centre`] work on sampled
//!   arrays supplied by a parser.
//!
//! Energies are in electronvolts and k-points are fractional reciprocal
//! coordinates unless a field documents otherwise.

mod bands;
mod dos;
mod mass;
mod validate;

pub use bands::{
    Band, BandGap, BandStructure, KPointLabel, Spin, SpinGap, band_gap, band_structure,
    path_distances, spin_band_structures, spin_gaps,
};
pub use dos::{
    DensityOfStates, DosBuffer, DosOptions, DosProjection, DosProjectionKind, Orbital,
    OrbitalWeight, Projection, d_band_centre, density_of_states, smear_dos, smear_levels,
};
pub use mass::{
    BandEdge, ELECTRON_MASS_FACTOR, EffectiveMass, SpinMasses, band_edge_masses, effective_mass,
    fit_parabolic_band, spin_band_edge_masses,
};

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

use validate::invalid;

/// Eigenvalues, occupations and projections of a calculation.
///
/// The first channel (`energies`, `occupations`) is spin up or unpolarized.
/// A spin-polarized calculation adds `spin_down` and, optionally, its
/// occupations. All bands of a channel share the same k-points, so
/// `energies[band][kpoint]` is the energy of one state.
///
/// Build it with [`ElectronicStructure::new`] and the `with_*` methods, then
/// call [`ElectronicStructure::validate`] once before handing it to the
/// analysis functions (they validate again, because the fields are public).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ElectronicStructure {
    /// Fractional k-point coordinates, one per column of `energies`.
    pub kpoints: Vec<[f64; 3]>,
    /// High-symmetry labels attached to k-point indices. Optional.
    #[serde(default)]
    pub labels: Vec<KPointLabel>,
    /// Band energies in eV, indexed as `energies[band][kpoint]`.
    pub energies: Vec<Vec<f64>>,
    /// Occupation numbers in `[0, 2]` with the same shape as `energies`.
    #[serde(default)]
    pub occupations: Option<Vec<Vec<f64>>>,
    /// Spin-down band energies in eV; same shape as `energies` when present.
    #[serde(default)]
    pub spin_down: Option<Vec<Vec<f64>>>,
    /// Spin-down occupations; requires `spin_down`.
    #[serde(default)]
    pub spin_down_occupations: Option<Vec<Vec<f64>>>,
    /// Projections of states onto elements and orbitals.
    #[serde(default)]
    pub projections: Vec<Projection>,
    /// Path coordinate for plotting and mass fits.
    ///
    /// Supply Cartesian reciprocal-space distances in inverse Angstrom to make
    /// effective masses physical. When absent, [`path_distances`] uses the
    /// Euclidean distance between the fractional k-points, which is only a
    /// plotting coordinate.
    #[serde(default)]
    pub distances: Option<Vec<f64>>,
    /// Fermi or reference energy in eV.
    ///
    /// Used to classify occupied states when `occupations` is absent.
    #[serde(default)]
    pub fermi: Option<f64>,
}

impl ElectronicStructure {
    /// Builds an unpolarized input from k-points and band energies.
    #[must_use]
    pub fn new(kpoints: Vec<[f64; 3]>, energies: Vec<Vec<f64>>) -> Self {
        Self {
            kpoints,
            labels: Vec::new(),
            energies,
            occupations: None,
            spin_down: None,
            spin_down_occupations: None,
            projections: Vec::new(),
            distances: None,
            fermi: None,
        }
    }

    /// Adds the occupations of the first channel.
    #[must_use]
    pub fn with_occupations(mut self, occupations: Vec<Vec<f64>>) -> Self {
        self.occupations = Some(occupations);
        self
    }

    /// Adds spin-down energies and enables the spin-resolved analysis.
    #[must_use]
    pub fn with_spin_down(mut self, energies: Vec<Vec<f64>>) -> Self {
        self.spin_down = Some(energies);
        self
    }

    /// Adds the spin-down occupations.
    #[must_use]
    pub fn with_spin_down_occupations(mut self, occupations: Vec<Vec<f64>>) -> Self {
        self.spin_down_occupations = Some(occupations);
        self
    }

    /// Adds element and orbital projections.
    #[must_use]
    pub fn with_projections(mut self, projections: Vec<Projection>) -> Self {
        self.projections = projections;
        self
    }

    /// Adds high-symmetry labels.
    #[must_use]
    pub fn with_labels(mut self, labels: Vec<KPointLabel>) -> Self {
        self.labels = labels;
        self
    }

    /// Sets the path coordinate used for plotting and mass fits.
    #[must_use]
    pub fn with_distances(mut self, distances: Vec<f64>) -> Self {
        self.distances = Some(distances);
        self
    }

    /// Sets the Fermi or reference energy in eV.
    #[must_use]
    pub fn with_fermi(mut self, fermi: f64) -> Self {
        self.fermi = Some(fermi);
        self
    }

    /// Number of k-points on the path.
    #[must_use]
    pub fn kpoint_count(&self) -> usize {
        self.kpoints.len()
    }

    /// Number of bands in the first channel.
    #[must_use]
    pub fn band_count(&self) -> usize {
        self.energies.len()
    }

    /// Number of spin channels: one when unpolarized, two when spin polarized.
    #[must_use]
    pub fn spin_count(&self) -> usize {
        1 + usize::from(self.spin_down.is_some())
    }

    /// Energies of one spin channel.
    #[must_use]
    pub fn energies_for(&self, spin: Spin) -> Option<&[Vec<f64>]> {
        match spin {
            Spin::Up => Some(&self.energies),
            Spin::Down => self.spin_down.as_deref(),
        }
    }

    /// Occupations of one spin channel, when reported by the calculation.
    #[must_use]
    pub fn occupations_for(&self, spin: Spin) -> Option<&[Vec<f64>]> {
        match spin {
            Spin::Up => self.occupations.as_deref(),
            Spin::Down => self.spin_down_occupations.as_deref(),
        }
    }

    /// Checks shapes, finiteness, label indices and projection indices.
    ///
    /// Returns [`Error::InvalidInput`] with a message that names the offending
    /// field. The analysis functions call this method before using the data.
    pub fn validate(&self) -> Result<(), Error> {
        if self.kpoints.is_empty() {
            return Err(invalid(
                "an electronic structure needs at least one k-point",
            ));
        }
        validate::kpoints(&self.kpoints)?;
        validate::energies(&self.energies, self.kpoints.len(), "energies")?;
        if let Some(occupations) = &self.occupations {
            validate::occupations(occupations, &self.energies, "occupations")?;
        }
        if let Some(spin_down) = &self.spin_down {
            if spin_down.len() != self.energies.len() {
                return Err(invalid(format!(
                    "spin_down has {} bands but energies has {}",
                    spin_down.len(),
                    self.energies.len()
                )));
            }
            validate::energies(spin_down, self.kpoints.len(), "spin_down")?;
            if let Some(occupations) = &self.spin_down_occupations {
                validate::occupations(occupations, spin_down, "spin_down_occupations")?;
            }
        } else if self.spin_down_occupations.is_some() {
            return Err(invalid(
                "spin_down_occupations requires spin_down band energies",
            ));
        }
        validate::labels(&self.labels, self.kpoints.len())?;
        if let Some(distances) = &self.distances {
            validate::distances(distances, self.kpoints.len(), "distances")?;
        }
        validate::projections(&self.projections, self.energies.len(), self.kpoints.len())?;
        validate::fermi(self.fermi)?;
        Ok(())
    }
}
