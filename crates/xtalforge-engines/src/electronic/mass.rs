//! Effective masses from parabolic fits around band edges.

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

use super::ElectronicStructure;
use super::bands::{BandStructure, Spin, spin_band_structures};
use super::validate::invalid;

/// `hbar^2 / m_e` in eV * Angstrom^2.
///
/// With the energy in eV and the wave-vector in inverse Angstrom, the
/// effective mass in electron masses is this constant divided by the second
/// derivative of the band.
pub const ELECTRON_MASS_FACTOR: f64 = 7.619_964_222_971_922;

/// Band edge at which a mass is evaluated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BandEdge {
    /// Conduction-band minimum; a normal band gives a positive mass.
    Conduction,
    /// Valence-band maximum; a normal band gives a negative (hole) mass.
    Valence,
}

impl BandEdge {
    /// Stable label used in reports and serialized output.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Conduction => "conduction",
            Self::Valence => "valence",
        }
    }
}

/// Effective mass and the parabola that produced it.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EffectiveMass {
    /// Signed mass in electron masses.
    ///
    /// Electron-like curvature gives a positive value; a valence-band maximum
    /// gives a negative value. Use [`EffectiveMass::magnitude`] for the
    /// absolute mass.
    pub mass: f64,
    /// Second derivative of the fitted band in eV * Angstrom^2.
    pub curvature: f64,
    /// Center of the fit in inverse Angstrom.
    pub k: f64,
    /// Energy at the center in eV.
    pub energy: f64,
    /// Number of k-points used by the least-squares fit.
    pub points: usize,
    /// Band edge that was fitted.
    pub edge: BandEdge,
}

impl EffectiveMass {
    /// Absolute mass in electron masses.
    #[must_use]
    pub fn magnitude(&self) -> f64 {
        self.mass.abs()
    }
}

/// Effective masses at the band edges of one spin channel.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SpinMasses {
    /// Spin channel.
    pub spin: Spin,
    /// Mass at the conduction-band minimum.
    pub electron: EffectiveMass,
    /// Mass at the valence-band maximum.
    pub hole: EffectiveMass,
}

/// Fits a parabola around one band edge and returns the mass in electron
/// masses.
///
/// `k` is the path coordinate in inverse Angstrom and must increase strictly.
/// The fit uses up to `window` neighbours on each side of `edge_index`; a
/// smaller one-sided window is used at the ends of the path, and at least three
/// points are required. The returned value is signed: a normal conduction band
/// gives a positive mass and a normal valence band a negative one.
pub fn effective_mass(
    k: &[f64],
    energies: &[f64],
    edge_index: usize,
    window: usize,
) -> Result<f64, Error> {
    Ok(fit_parabolic_band(k, energies, edge_index, window, BandEdge::Conduction)?.mass)
}

/// Fits a parabola around one band edge and returns the full result.
///
/// See [`effective_mass`] for the units and the window rules.
pub fn fit_parabolic_band(
    k: &[f64],
    energies: &[f64],
    edge_index: usize,
    window: usize,
    edge: BandEdge,
) -> Result<EffectiveMass, Error> {
    if k.len() != energies.len() {
        return Err(invalid(format!(
            "the fit has {} k-points but {} energies",
            k.len(),
            energies.len()
        )));
    }
    if k.len() < 3 {
        return Err(invalid("a parabolic fit needs at least three k-points"));
    }
    if window == 0 {
        return Err(invalid("the fit window must be at least one k-point"));
    }
    if edge_index >= k.len() {
        return Err(invalid(format!(
            "edge index {edge_index} is out of range for {} k-points",
            k.len()
        )));
    }
    if k.iter().any(|value| !value.is_finite()) {
        return Err(invalid("the k-point coordinates are not finite"));
    }
    if energies.iter().any(|value| !value.is_finite()) {
        return Err(invalid("the band energies are not finite"));
    }
    for pair in k.windows(2) {
        if let [left, right] = pair
            && right <= left
        {
            return Err(invalid("the k-point coordinates must increase strictly"));
        }
    }
    let start = edge_index.saturating_sub(window);
    let end = edge_index
        .saturating_add(window)
        .saturating_add(1)
        .min(k.len());
    let count = end - start;
    if count < 3 {
        return Err(invalid("the fit window contains fewer than three k-points"));
    }
    let center = *k
        .get(edge_index)
        .ok_or_else(|| invalid("the edge index is out of range"))?;
    let energy = *energies
        .get(edge_index)
        .ok_or_else(|| invalid("the edge index is out of range"))?;

    let mut count_f = 0.0;
    let mut s1 = 0.0;
    let mut s2 = 0.0;
    let mut s3 = 0.0;
    let mut s4 = 0.0;
    let mut sy = 0.0;
    let mut sxy = 0.0;
    let mut sx2y = 0.0;
    for (x, y) in k
        .iter()
        .skip(start)
        .take(count)
        .zip(energies.iter().skip(start))
    {
        let dx = x - center;
        let dx2 = dx * dx;
        count_f += 1.0;
        s1 += dx;
        s2 += dx2;
        s3 += dx2 * dx;
        s4 += dx2 * dx2;
        sy += y;
        sxy += dx * y;
        sx2y += dx2 * y;
    }

    // Normal equations of `E = a + b dx + c dx^2`, solved with Cramer's rule.
    let determinant =
        count_f * (s2 * s4 - s3 * s3) - s1 * (s1 * s4 - s3 * s2) + s2 * (s1 * s3 - s2 * s2);
    let scale = count_f * s2 * s4;
    if !determinant.is_finite()
        || !scale.is_finite()
        || scale <= 0.0
        || determinant.abs() < 1e-14 * scale
    {
        return Err(invalid(
            "the k-points are too clustered for a parabolic fit",
        ));
    }
    let determinant_c =
        count_f * (s2 * sx2y - s3 * sxy) - s1 * (s1 * sx2y - s3 * sy) + s2 * (s1 * sxy - s2 * sy);
    let coefficient = determinant_c / determinant;
    let curvature = 2.0 * coefficient;
    let mass = ELECTRON_MASS_FACTOR / curvature;
    if !curvature.is_finite() || !mass.is_finite() {
        return Err(invalid(
            "the fitted band is flat and has no finite effective mass",
        ));
    }
    Ok(EffectiveMass {
        mass,
        curvature,
        k: center,
        energy,
        points: count,
        edge,
    })
}

/// Effective masses at the valence and conduction edges of one channel.
///
/// The returned vector holds the electron mass first and the hole mass second.
/// The edge classification needs occupations or a Fermi energy; see
/// [`BandStructure::gap`]. Distances must be Cartesian reciprocal-space
/// distances in inverse Angstrom; see
/// [`ElectronicStructure::distances`](super::ElectronicStructure::distances).
pub fn band_edge_masses(
    structure: &BandStructure,
    window: usize,
) -> Result<Vec<EffectiveMass>, Error> {
    let (electron, hole) = edge_masses(structure, window)?;
    Ok(vec![electron, hole])
}

/// Effective masses at the band edges of every spin channel.
pub fn spin_band_edge_masses(
    input: &ElectronicStructure,
    window: usize,
) -> Result<Vec<SpinMasses>, Error> {
    spin_band_structures(input)?
        .into_iter()
        .map(|structure| {
            let spin = structure.spin;
            let (electron, hole) = edge_masses(&structure, window)?;
            Ok(SpinMasses {
                spin,
                electron,
                hole,
            })
        })
        .collect()
}

/// Fits the states at the valence and conduction edges of one channel.
fn edge_masses(
    structure: &BandStructure,
    window: usize,
) -> Result<(EffectiveMass, EffectiveMass), Error> {
    let gap = structure.gap()?;
    let electron_band = edge_band(structure, gap.cbm_kpoint, gap.cbm)?;
    let hole_band = edge_band(structure, gap.vbm_kpoint, gap.vbm)?;
    let electron = fit_selected(
        structure,
        electron_band,
        gap.cbm_kpoint,
        window,
        BandEdge::Conduction,
    )?;
    let hole = fit_selected(
        structure,
        hole_band,
        gap.vbm_kpoint,
        window,
        BandEdge::Valence,
    )?;
    Ok((electron, hole))
}

/// Finds the band that carries the requested edge energy at a k-point.
fn edge_band(structure: &BandStructure, kpoint: usize, target: f64) -> Result<usize, Error> {
    let mut best: Option<(usize, f64)> = None;
    for (index, band) in structure.bands.iter().enumerate() {
        let energy = band
            .energies
            .get(kpoint)
            .copied()
            .ok_or_else(|| invalid(format!("band {index} is missing k-point {kpoint}")))?;
        let distance = (energy - target).abs();
        if best.is_none_or(|(_, current)| distance < current) {
            best = Some((index, distance));
        }
    }
    best.map(|(index, _)| index)
        .ok_or_else(|| invalid("the band structure has no bands"))
}

/// Fits one selected band at a k-point.
fn fit_selected(
    structure: &BandStructure,
    band: usize,
    kpoint: usize,
    window: usize,
    edge: BandEdge,
) -> Result<EffectiveMass, Error> {
    let energies = structure
        .bands
        .get(band)
        .map(|band| band.energies.as_slice())
        .ok_or_else(|| invalid(format!("band {band} is out of range")))?;
    fit_parabolic_band(&structure.distances, energies, kpoint, window, edge)
}
