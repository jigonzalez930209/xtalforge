//! Phonon density of states and inelastic-neutron weights (Task 19.4.3).
//!
//! [`phonon_dos`] rebuilds the Gaussian-smeared DOS of a [`PhononResult`]
//! without requiring `PhononOptions::dos`, using the same histogram policy as
//! the q-mesh code: bins between the data minimum and maximum (or the caller's
//! range), unit-height Gaussian smearing and unit-area normalisation.
//!
//! [`inelastic_neutron_dos`] weights each mode by the coherent inelastic
//! scattering weight
//!
//! ```text
//! w_m = sum_i (b_i^2 / m_i) |u_{i,m}|^2 / sum_i |u_{i,m}|^2
//! ```
//!
//! in Angstrom^2 / amu, where `b_i` is the bound coherent neutron scattering
//! length, `m_i` the atomic mass and `u_{i,m}` the Cartesian displacement of
//! atom `i` in mode `m`.  The ratio removes the arbitrary mode normalisation
//! used by [`VibrationalModes`], so `w_m` is unchanged when a mode is scaled.
//! The weighted DOS is normalised to unit area like the plain one; the
//! absolute scattering scale (instrument factors, Debye-Waller attenuation and
//! the energy-dependent detector response) is not modelled.  The per-branch
//! weights are the Gamma-point eigenvector weights; a q-point frequency is
//! assigned to the branch with the same index in the ascending eigenvalue list,
//! so a branch crossing mixes Gamma weights.  That approximation is documented
//! rather than hidden and is exact when branches do not cross or when the
//! weights are uniform.
//!
//! **Dependency choice.**  `xtalforge-ff` deliberately does not depend on
//! `xtalforge-analysis`; doing so would pull the I/O and symmetry crates and
//! `moyo` into the force-field layer.  [`NeutronWeights`] therefore accepts the
//! scattering lengths as a table, a callback or an explicit vector, and a
//! caller in the facade (which links both crates) can fill it from
//! `xtalforge_analysis::scattering::neutron_scattering_length`.  Negative
//! lengths (H, Ti, Mn, V) are kept because the weight uses their square.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Structure, element_info};

use crate::Error;

use super::modes::VibrationalModes;
use super::phonons::PhononResult;
use super::qmesh::{DosOptions, PhononDos};

/// Bound coherent neutron scattering lengths of one P1 cell, in Angstrom.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NeutronWeights {
    /// One length per P1 atom, in Angstrom; negative values are physical.
    pub lengths: Vec<f64>,
}

impl NeutronWeights {
    /// Wraps explicit lengths.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when a length is not finite.
    pub fn from_lengths(lengths: Vec<f64>) -> Result<Self, Error> {
        if lengths.iter().any(|length| !length.is_finite()) {
            return Err(Error::InvalidInput(
                "neutron scattering lengths must be finite".to_owned(),
            ));
        }
        Ok(Self { lengths })
    }

    /// Looks the lengths up in a caller-supplied table keyed by element symbol.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when an element is missing from the
    /// table, and [`Error::UnknownElement`] for an unknown site element.
    pub fn from_table(structure: &Structure, table: &BTreeMap<String, f64>) -> Result<Self, Error> {
        Self::from_callback(structure, |symbol| {
            table.get(symbol).copied().ok_or_else(|| {
                Error::InvalidInput(format!("no neutron scattering length for element {symbol}"))
            })
        })
    }

    /// Builds the lengths with a callback that receives each element symbol.
    ///
    /// This is the bridge to an external table such as
    /// `xtalforge_analysis::scattering::neutron_scattering_length` without
    /// making this crate depend on that one.
    ///
    /// # Errors
    ///
    /// Returns the callback error for an element without a length,
    /// [`Error::UnknownElement`] for an unknown site element, and
    /// [`Error::InvalidInput`] for a non-finite result.
    pub fn from_callback<F>(structure: &Structure, mut callback: F) -> Result<Self, Error>
    where
        F: FnMut(&str) -> Result<f64, Error>,
    {
        let mut lengths = Vec::with_capacity(structure.sites().len());
        for site in structure.sites() {
            let element = element_info(site.element())?;
            let length = callback(&element.symbol)?;
            if !length.is_finite() {
                return Err(Error::InvalidInput(format!(
                    "neutron scattering length for {} is not finite",
                    element.symbol
                )));
            }
            lengths.push(length);
        }
        Ok(Self { lengths })
    }

    /// Number of stored lengths.
    #[must_use]
    pub fn len(&self) -> usize {
        self.lengths.len()
    }

    /// Whether no length is stored.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.lengths.is_empty()
    }
}

/// Plain and neutron-weighted phonon DOS of one [`PhononResult`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InelasticNeutronDos {
    /// Unweighted DOS, normalised to unit area.
    pub dos: PhononDos,
    /// `b^2 / m`-weighted DOS, normalised to unit area.
    pub weighted: PhononDos,
    /// Gamma-point weight of each mode in Angstrom^2 / amu.
    pub mode_weights: Vec<f64>,
    /// Sum of `mode_weights`.
    pub total_weight: f64,
}

/// Rebuilds the phonon DOS of a result, whether or not it carries one.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid DOS options, an empty frequency
/// list or a non-finite frequency.
pub fn phonon_dos(result: &PhononResult, options: &DosOptions) -> Result<PhononDos, Error> {
    let samples = result
        .frequencies
        .iter()
        .flatten()
        .copied()
        .map(|frequency| (frequency, 1.0))
        .collect::<Vec<_>>();
    histogram(&samples, options)
}

/// Gamma-point neutron weight of every mode in a [`VibrationalModes`].
///
/// The result has one entry per mode (including zero and imaginary modes) and
/// is in Angstrom^2 / amu.  A mode whose Cartesian norm is zero gets weight
/// zero instead of a division by zero.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-P1 structure, an atom-count
/// mismatch, a malformed mode list, a non-positive mass, or a non-finite
/// displacement.
pub fn mode_neutron_weights(
    structure: &Structure,
    modes: &VibrationalModes,
    weights: &NeutronWeights,
) -> Result<Vec<f64>, Error> {
    if !structure.is_p1() {
        return Err(Error::InvalidInput(
            "neutron weights need a P1 structure; expand the asymmetric unit first".to_owned(),
        ));
    }
    let count = structure.sites().len();
    if modes.n_atoms() != count {
        return Err(Error::InvalidInput(format!(
            "modes describe {} atoms but the structure has {count}",
            modes.n_atoms()
        )));
    }
    if weights.lengths.len() != count {
        return Err(Error::InvalidInput(format!(
            "got {} neutron lengths for {count} atoms",
            weights.lengths.len()
        )));
    }
    if modes.frequencies.len() != modes.modes.len() {
        return Err(Error::InvalidInput(format!(
            "modes store {} frequencies but {} displacement fields",
            modes.frequencies.len(),
            modes.modes.len()
        )));
    }
    let mut masses = Vec::with_capacity(count);
    for site in structure.sites() {
        let element = element_info(site.element())?;
        if element.mass <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "element {} has no positive mass",
                site.element()
            )));
        }
        masses.push(element.mass);
    }

    let mut mode_weights = Vec::with_capacity(modes.n_modes());
    for (mode, displacements) in modes.modes.iter().enumerate() {
        if displacements.len() != count {
            return Err(Error::InvalidInput(format!(
                "mode {mode} has {} atoms, expected {count}",
                displacements.len()
            )));
        }
        let mut norm = 0.0;
        let mut weighted = 0.0;
        for (index, point) in displacements.iter().enumerate() {
            if point.iter().any(|value| !value.is_finite()) {
                return Err(Error::InvalidInput(format!(
                    "mode {mode} has a non-finite displacement"
                )));
            }
            let length_squared = point[0] * point[0] + point[1] * point[1] + point[2] * point[2];
            norm += length_squared;
            let Some(length) = weights.lengths.get(index).copied() else {
                continue;
            };
            let Some(mass) = masses.get(index).copied() else {
                continue;
            };
            weighted += length * length / mass * length_squared;
        }
        mode_weights.push(if norm > 0.0 { weighted / norm } else { 0.0 });
    }
    Ok(mode_weights)
}

/// Builds the plain and neutron-weighted DOS of a phonon result.
///
/// `structure` must be the P1 cell the result was calculated for, because the
/// per-atom scattering lengths and masses are matched by atom index.
///
/// # Errors
///
/// Returns the errors of [`mode_neutron_weights`] and [`phonon_dos`], and
/// [`Error::InvalidInput`] when a q-point frequency list does not match the
/// number of Gamma modes or when the weights sum to zero.
pub fn inelastic_neutron_dos(
    structure: &Structure,
    result: &PhononResult,
    weights: &NeutronWeights,
    options: &DosOptions,
) -> Result<InelasticNeutronDos, Error> {
    let mode_weights = mode_neutron_weights(structure, &result.modes, weights)?;
    let mut samples = Vec::new();
    for (q_index, values) in result.frequencies.iter().enumerate() {
        if values.len() != mode_weights.len() {
            return Err(Error::InvalidInput(format!(
                "q-point {q_index} has {} frequencies but the Gamma modes have {}",
                values.len(),
                mode_weights.len()
            )));
        }
        for (mode, value) in values.iter().enumerate() {
            let weight = mode_weights.get(mode).copied().unwrap_or(0.0);
            samples.push((*value, weight));
        }
    }
    let dos = phonon_dos(result, options)?;
    let weighted = histogram(&samples, options)?;
    let total_weight = mode_weights.iter().sum();
    Ok(InelasticNeutronDos {
        dos,
        weighted,
        mode_weights,
        total_weight,
    })
}

fn histogram(samples: &[(f64, f64)], options: &DosOptions) -> Result<PhononDos, Error> {
    options.validate()?;
    if samples.is_empty() {
        return Err(Error::InvalidInput(
            "phonon DOS needs at least one frequency".to_owned(),
        ));
    }
    if samples
        .iter()
        .any(|(frequency, weight)| !frequency.is_finite() || !weight.is_finite())
    {
        return Err(Error::InvalidInput(
            "phonon DOS frequencies and weights must be finite".to_owned(),
        ));
    }
    let first = samples
        .iter()
        .map(|(frequency, _)| *frequency)
        .fold(f64::INFINITY, f64::min);
    let last = samples
        .iter()
        .map(|(frequency, _)| *frequency)
        .fold(f64::NEG_INFINITY, f64::max);
    let min = options.min.unwrap_or(first);
    let max = options.max.unwrap_or(last.max(min + 1.0));
    if min >= max {
        return Err(Error::InvalidInput("phonon DOS range is empty".to_owned()));
    }
    let bins = options.bins;
    let step = (max - min) / bins as f64;
    let mut density = vec![0.0; bins];
    let sigma2 = 2.0 * options.sigma * options.sigma;
    for (frequency, weight) in samples {
        for (index, cell) in density.iter_mut().enumerate() {
            let center = min + (index as f64 + 0.5) * step;
            let offset = center - frequency;
            *cell += weight * (-offset * offset / sigma2).exp();
        }
    }
    let integral = density.iter().sum::<f64>() * step;
    if integral <= 0.0 {
        return Err(Error::InvalidInput(
            "phonon DOS has zero total weight".to_owned(),
        ));
    }
    for value in &mut density {
        *value /= integral;
    }
    Ok(PhononDos {
        frequencies: (0..bins)
            .map(|index| min + (index as f64 + 0.5) * step)
            .collect(),
        density,
        sigma: options.sigma,
    })
}
