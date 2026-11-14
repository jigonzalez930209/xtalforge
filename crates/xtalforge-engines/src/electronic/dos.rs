//! Total and projected densities of states.

mod smear;

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

use super::ElectronicStructure;
use super::validate::invalid;
use smear::{integrate, validate_sigma};

pub use smear::{d_band_centre, smear_dos, smear_levels};

/// Square root of two pi, used by the Gaussian kernel.
const SQRT_TWO_PI: f64 = 2.506_628_274_631_000_7;

/// Orbital character of a projection.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Orbital {
    /// s orbital.
    S,
    /// p orbital.
    P,
    /// d orbital, the one used by the d-band centre.
    D,
    /// f orbital.
    F,
    /// Any other label, such as `dxy` or `pz`.
    Other(String),
}

impl Orbital {
    /// Stable label used in reports and serialized output.
    #[must_use]
    pub fn as_str(&self) -> &str {
        match self {
            Self::S => "s",
            Self::P => "p",
            Self::D => "d",
            Self::F => "f",
            Self::Other(label) => label,
        }
    }

    /// Parses an orbital label, keeping unknown labels verbatim.
    #[must_use]
    pub fn parse(label: &str) -> Self {
        let trimmed = label.trim();
        match trimmed.to_ascii_lowercase().as_str() {
            "s" => Self::S,
            "p" => Self::P,
            "d" => Self::D,
            "f" => Self::F,
            _ => Self::Other(trimmed.to_owned()),
        }
    }
}

/// One element/orbital weight of a projected state.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OrbitalWeight {
    /// Site index in the originating structure.
    pub site: usize,
    /// Element symbol.
    pub element: String,
    /// Orbital character.
    pub orbital: Orbital,
    /// Projection weight of the state; finite and non-negative.
    pub weight: f64,
}

/// Projection of one band state onto element/orbital weights.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Projection {
    /// Band index.
    pub band: usize,
    /// K-point index.
    pub kpoint: usize,
    /// Weights of this state. An empty list contributes nothing.
    #[serde(default)]
    pub weights: Vec<OrbitalWeight>,
}

/// Kind of a DOS projection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DosProjectionKind {
    /// Sum over all sites of one element.
    Element,
    /// Sum over all sites with one orbital character.
    Orbital,
}

/// One projected DOS channel on the energy grid.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DosProjection {
    /// Label shown in a plot: an element symbol or an orbital name.
    pub label: String,
    /// Kind of grouping.
    pub kind: DosProjectionKind,
    /// Element symbol for an element projection.
    #[serde(default)]
    pub element: Option<String>,
    /// Orbital for an orbital projection.
    #[serde(default)]
    pub orbital: Option<Orbital>,
    /// Values on the energy grid of the parent [`DensityOfStates`].
    pub values: Vec<f64>,
}

impl DosProjection {
    /// Typed buffer for plotting on the bindings.
    #[must_use]
    pub fn buffer(&self, energies: &[f64]) -> DosBuffer {
        DosBuffer {
            energies: energies.to_vec(),
            values: self.values.clone(),
        }
    }
}

/// Energy grid and values of one DOS channel.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DosBuffer {
    /// Energy grid in eV.
    pub energies: Vec<f64>,
    /// DOS values on the grid.
    pub values: Vec<f64>,
}

/// Options for [`density_of_states`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DosOptions {
    /// Standard deviation of the Gaussian in eV; must be positive.
    pub sigma: f64,
    /// Number of points on the energy grid; at least two.
    pub points: usize,
    /// Explicit energy range in eV. When absent the data range is padded by
    /// five standard deviations on both sides.
    #[serde(default)]
    pub energy_range: Option<(f64, f64)>,
}

impl Default for DosOptions {
    fn default() -> Self {
        Self {
            sigma: 0.1,
            points: 1001,
            energy_range: None,
        }
    }
}

impl DosOptions {
    /// Options with a Gaussian width and the default grid.
    #[must_use]
    pub fn new(sigma: f64) -> Self {
        Self {
            sigma,
            ..Self::default()
        }
    }

    /// Checks the width, the grid size and the explicit range.
    pub fn validate(&self) -> Result<(), Error> {
        validate_sigma(self.sigma)?;
        if self.points < 2 {
            return Err(invalid("a DOS grid needs at least two points"));
        }
        if let Some((min, max)) = self.energy_range {
            if !min.is_finite() || !max.is_finite() {
                return Err(invalid("the DOS energy range is not finite"));
            }
            if min >= max {
                return Err(invalid(format!(
                    "the DOS energy range must increase, got [{min}, {max}]"
                )));
            }
        }
        Ok(())
    }
}

/// Total and projected density of states on an energy grid.
///
/// Values are a sum over the supplied k-points and bands without k-point
/// weights: the integral is the number of states when the input is a uniform
/// mesh and the caller accounts for the mesh normalization.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DensityOfStates {
    /// Energy grid in eV.
    pub energies: Vec<f64>,
    /// Total DOS, both spin channels when present.
    pub total: Vec<f64>,
    /// Element projections (sorted by symbol) followed by orbital projections
    /// (sorted by orbital).
    pub projections: Vec<DosProjection>,
    /// Fermi or reference energy copied from the input.
    #[serde(default)]
    pub fermi: Option<f64>,
}

impl DensityOfStates {
    /// Typed buffer of the total DOS.
    #[must_use]
    pub fn total_buffer(&self) -> DosBuffer {
        DosBuffer {
            energies: self.energies.clone(),
            values: self.total.clone(),
        }
    }

    /// First projection with a label, whatever its kind.
    ///
    /// An element and an orbital can share a label; use [`Self::element`] or
    /// [`Self::orbital`] to select a specific kind.
    #[must_use]
    pub fn projection(&self, label: &str) -> Option<&DosProjection> {
        self.projections
            .iter()
            .find(|projection| projection.label == label)
    }

    /// Element projection of a symbol.
    #[must_use]
    pub fn element(&self, element: &str) -> Option<&DosProjection> {
        self.projections.iter().find(|projection| {
            projection.kind == DosProjectionKind::Element
                && projection.element.as_deref() == Some(element)
        })
    }

    /// Orbital projection of one orbital character.
    #[must_use]
    pub fn orbital(&self, orbital: Orbital) -> Option<&DosProjection> {
        self.projections.iter().find(|projection| {
            projection.kind == DosProjectionKind::Orbital
                && projection.orbital.as_ref() == Some(&orbital)
        })
    }

    /// Integral of the total DOS over the grid, in states.
    pub fn total_states(&self) -> Result<f64, Error> {
        integrate(&self.energies, &self.total)
    }

    /// Centre of the summed d-orbital projection in eV.
    pub fn d_band_centre(&self) -> Result<f64, Error> {
        let mut values = vec![0.0; self.energies.len()];
        let mut found = false;
        for projection in self.projections.iter().filter(|projection| {
            projection.kind == DosProjectionKind::Orbital
                && projection.orbital.as_ref() == Some(&Orbital::D)
        }) {
            found = true;
            for (slot, value) in values.iter_mut().zip(&projection.values) {
                *slot += value;
            }
        }
        if !found {
            return Err(invalid("the DOS has no d-orbital projection"));
        }
        d_band_centre(&self.energies, &values)
    }
}

/// Builds total, element- and orbital-projected DOS from an electronic
/// structure.
///
/// The total channel contains every band of every spin channel. Projections
/// add the Gaussian of each projected state weighted by its element/orbital
/// weight, so a weight of one contributes a full state to that channel.
pub fn density_of_states(
    input: &ElectronicStructure,
    options: &DosOptions,
) -> Result<DensityOfStates, Error> {
    input.validate()?;
    options.validate()?;
    let mut levels = Vec::new();
    for row in &input.energies {
        levels.extend(row.iter().copied());
    }
    if let Some(spin_down) = &input.spin_down {
        for row in spin_down {
            levels.extend(row.iter().copied());
        }
    }
    if levels.is_empty() {
        return Err(invalid("the electronic structure has no energies"));
    }
    let (min, max) = match options.energy_range {
        Some(range) => range,
        None => {
            let mut min = f64::INFINITY;
            let mut max = f64::NEG_INFINITY;
            for level in &levels {
                min = min.min(*level);
                max = max.max(*level);
            }
            let padding = 5.0 * options.sigma;
            (min - padding, max + padding)
        }
    };
    let points = options.points;
    let span = max - min;
    let step = span / (points - 1) as f64;
    let energies: Vec<f64> = (0..points).map(|index| min + step * index as f64).collect();
    let total = smear_levels(&levels, &energies, options.sigma)?;
    let projections = projected_dos(input, &energies, options.sigma)?;
    Ok(DensityOfStates {
        energies,
        total,
        projections,
        fermi: input.fermi,
    })
}

/// Builds element and orbital projections on an existing grid.
fn projected_dos(
    input: &ElectronicStructure,
    grid: &[f64],
    sigma: f64,
) -> Result<Vec<DosProjection>, Error> {
    let mut elements: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    let mut orbitals: BTreeMap<Orbital, Vec<f64>> = BTreeMap::new();
    for projection in &input.projections {
        let level = input
            .energies
            .get(projection.band)
            .and_then(|row| row.get(projection.kpoint))
            .copied()
            .ok_or_else(|| invalid("a projection refers to a missing state"))?;
        let values = gaussian_values(level, grid, sigma);
        for weight in &projection.weights {
            let element = elements
                .entry(weight.element.clone())
                .or_insert_with(|| vec![0.0; grid.len()]);
            add_scaled(element, &values, weight.weight);
            let orbital = orbitals
                .entry(weight.orbital.clone())
                .or_insert_with(|| vec![0.0; grid.len()]);
            add_scaled(orbital, &values, weight.weight);
        }
    }
    let mut projections = Vec::with_capacity(elements.len() + orbitals.len());
    for (element, values) in elements {
        projections.push(DosProjection {
            label: element.clone(),
            kind: DosProjectionKind::Element,
            element: Some(element),
            orbital: None,
            values,
        });
    }
    for (orbital, values) in orbitals {
        projections.push(DosProjection {
            label: orbital.as_str().to_owned(),
            kind: DosProjectionKind::Orbital,
            element: None,
            orbital: Some(orbital),
            values,
        });
    }
    Ok(projections)
}

/// Adds a scaled channel to an accumulator.
fn add_scaled(target: &mut [f64], source: &[f64], scale: f64) {
    for (slot, value) in target.iter_mut().zip(source) {
        *slot += value * scale;
    }
}

/// Normalized Gaussian of `center` sampled on `grid`.
fn gaussian_values(center: f64, grid: &[f64], sigma: f64) -> Vec<f64> {
    let normalization = 1.0 / (sigma * SQRT_TWO_PI);
    grid.iter()
        .map(|energy| {
            let scaled = (energy - center) / sigma;
            normalization * (-0.5 * scaled * scaled).exp()
        })
        .collect()
}
