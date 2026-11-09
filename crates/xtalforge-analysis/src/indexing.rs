//! Powder indexing: metric families and figures of merit (Task 7.5).
//!
//! [`index_cubic`] is the original single-parameter search. [`index_pattern`]
//! generalises it to the tetragonal, hexagonal, orthorhombic and monoclinic
//! metric families, ranks the candidates with the de Wolff M20 and
//! Smith-Snyder F_N figures of merit, and resolves a cell found by several
//! families in favour of the highest symmetry.

mod cubic;
mod fom;
mod general;

pub use cubic::{CubicIndex, index_cubic};
pub use fom::{
    FiguresOfMerit, figures_of_merit, match_peaks, peak_q, predict_reflections, two_theta_of_q,
};

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Lattice};
use xtalforge_symmetry::{Reduction, reduce_lattice};

use crate::peaks::Peak;

/// One indexed reflection.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct IndexedPeak {
    /// Miller indices.
    pub hkl: [i32; 3],
    /// Peak position in degrees two-theta.
    pub two_theta: f64,
}

/// Metric family of a cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LatticeSystem {
    /// `a = b = c`, all angles 90 degrees.
    Cubic,
    /// `a = b`, all angles 90 degrees.
    Tetragonal,
    /// `a = b`, `gamma = 120` degrees.
    Hexagonal,
    /// Three lengths, all angles 90 degrees.
    Orthorhombic,
    /// Unique `b` axis, `beta` free.
    Monoclinic,
}

impl LatticeSystem {
    /// Stable name.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::Cubic => "cubic",
            Self::Tetragonal => "tetragonal",
            Self::Hexagonal => "hexagonal",
            Self::Orthorhombic => "orthorhombic",
            Self::Monoclinic => "monoclinic",
        }
    }

    /// Number of metric parameters the family leaves free.
    #[must_use]
    pub const fn parameters(self) -> usize {
        match self {
            Self::Cubic => 1,
            Self::Tetragonal | Self::Hexagonal => 2,
            Self::Orthorhombic => 3,
            Self::Monoclinic => 4,
        }
    }

    /// Builds a cell from the solved metric parameters.
    #[must_use]
    pub fn cell(self, values: &[f64]) -> Option<CellMetric> {
        let value = |index: usize| values.get(index).copied();
        match self {
            Self::Cubic => {
                let a = 1.0 / value(0)?.max(0.0).sqrt();
                Some(CellMetric::new(a, a, a, 90.0, 90.0, 90.0))
            }
            Self::Tetragonal => {
                let a = 1.0 / value(0)?.max(0.0).sqrt();
                let c = 1.0 / value(1)?.max(0.0).sqrt();
                Some(CellMetric::new(a, a, c, 90.0, 90.0, 90.0))
            }
            Self::Hexagonal => {
                let a = (4.0 / (3.0 * value(0)?.max(0.0))).sqrt();
                let c = 1.0 / value(1)?.max(0.0).sqrt();
                Some(CellMetric::new(a, a, c, 90.0, 90.0, 120.0))
            }
            Self::Orthorhombic => {
                let a = 1.0 / value(0)?.max(0.0).sqrt();
                let b = 1.0 / value(1)?.max(0.0).sqrt();
                let c = 1.0 / value(2)?.max(0.0).sqrt();
                Some(CellMetric::new(a, b, c, 90.0, 90.0, 90.0))
            }
            Self::Monoclinic => {
                let a_metric = value(0)?;
                let b_metric = value(1)?;
                let c_metric = value(2)?;
                let cross = value(3)?;
                if a_metric <= 0.0 || b_metric <= 0.0 || c_metric <= 0.0 {
                    return None;
                }
                let cosine = -cross / (2.0 * (a_metric * c_metric).sqrt());
                if cosine.abs() >= 1.0 {
                    return None;
                }
                let beta = cosine.acos().to_degrees();
                let sine = cosine.acos().sin();
                if sine.abs() < 1e-9 {
                    return None;
                }
                let a = 1.0 / (a_metric * sine * sine).sqrt();
                let b = 1.0 / b_metric.sqrt();
                let c = 1.0 / (c_metric * sine * sine).sqrt();
                Some(CellMetric::new(a, b, c, 90.0, beta, 90.0))
            }
        }
    }
}

/// Cell parameters in Angstrom and degrees.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CellMetric {
    /// `a` length.
    pub a: f64,
    /// `b` length.
    pub b: f64,
    /// `c` length.
    pub c: f64,
    /// Angle between `b` and `c`.
    pub alpha: f64,
    /// Angle between `a` and `c`.
    pub beta: f64,
    /// Angle between `a` and `b`.
    pub gamma: f64,
}

impl CellMetric {
    /// Creates a cell; no validation is applied.
    #[must_use]
    pub const fn new(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> Self {
        Self {
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
        }
    }

    /// Unit-cell volume in cubic Angstrom.
    #[must_use]
    pub fn volume(&self) -> f64 {
        let alpha = self.alpha.to_radians();
        let beta = self.beta.to_radians();
        let gamma = self.gamma.to_radians();
        let cosine = 1.0 - alpha.cos().powi(2) - beta.cos().powi(2) - gamma.cos().powi(2)
            + 2.0 * alpha.cos() * beta.cos() * gamma.cos();
        self.a * self.b * self.c * cosine.max(0.0).sqrt()
    }

    /// `Q = 1/d^2` of a reflection, or `None` for the 000 reflection.
    ///
    /// Uses the reciprocal metric tensor, so every family (including the
    /// oblique ones) takes the same path.
    #[must_use]
    pub fn reciprocal_squared(&self, hkl: [i32; 3]) -> Option<f64> {
        if hkl == [0, 0, 0] {
            return None;
        }
        let [h, k, l] = hkl.map(f64::from);
        let alpha = self.alpha.to_radians();
        let beta = self.beta.to_radians();
        let gamma = self.gamma.to_radians();
        let volume = self.volume();
        if volume <= 0.0 {
            return None;
        }
        let s11 = self.b * self.b * self.c * self.c * alpha.sin().powi(2);
        let s22 = self.a * self.a * self.c * self.c * beta.sin().powi(2);
        let s33 = self.a * self.a * self.b * self.b * gamma.sin().powi(2);
        let s12 = self.a * self.b * self.c * self.c * (alpha.cos() * beta.cos() - gamma.cos());
        let s13 = self.a * self.b * self.b * self.c * (alpha.cos() * gamma.cos() - beta.cos());
        let s23 = self.a * self.a * self.b * self.c * (beta.cos() * gamma.cos() - alpha.cos());
        let q = (h * h * s11
            + k * k * s22
            + l * l * s33
            + 2.0 * h * k * s12
            + 2.0 * h * l * s13
            + 2.0 * k * l * s23)
            / (volume * volume);
        (q > 0.0 && q.is_finite()).then_some(q)
    }

    /// Rounded key used to merge the same cell found by several families.
    #[must_use]
    pub fn key(&self) -> [i64; 6] {
        [
            (self.a * 1000.0).round() as i64,
            (self.b * 1000.0).round() as i64,
            (self.c * 1000.0).round() as i64,
            (self.alpha * 100.0).round() as i64,
            (self.beta * 100.0).round() as i64,
            (self.gamma * 100.0).round() as i64,
        ]
    }

    /// Whether every parameter is finite, inside the search bounds and the
    /// volume is not absurd.
    #[must_use]
    pub fn is_plausible(&self, max_volume: f64) -> bool {
        let lengths = [self.a, self.b, self.c];
        let angles = [self.alpha, self.beta, self.gamma];
        lengths.iter().all(|value| (2.0..=60.0).contains(value))
            && angles.iter().all(|value| (30.0..=150.0).contains(value))
            && self.volume() >= 10.0
            && self.volume() <= max_volume
    }
}

/// One ranked indexing candidate.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IndexedCell {
    /// Metric family.
    pub system: LatticeSystem,
    /// Refined cell parameters.
    pub cell: CellMetric,
    /// Indexed peaks.
    pub peaks: Vec<IndexedPeak>,
    /// Fraction of observed peaks that found an integer hkl.
    pub completeness: f64,
    /// Mean absolute observed minus calculated 2-theta in degrees.
    pub mean_error: f64,
    /// de Wolff M20 figure of merit (higher is better).
    pub m20: f64,
    /// Smith-Snyder F_N figure of merit (higher is better).
    pub f_n: f64,
    /// Fraction of the predicted lines below `Q20` that are observed.
    pub coverage: f64,
}

/// Search settings for [`index_pattern`].
#[derive(Debug, Clone, PartialEq)]
pub struct IndexingOptions {
    /// Metric families to try.
    pub systems: Vec<LatticeSystem>,
    /// Largest |h|, |k| or |l| enumerated for the solving reflections.
    pub max_index: i32,
    /// Position tolerance in degrees two-theta.
    pub tolerance: f64,
    /// Largest accepted cell volume in cubic Angstrom.
    pub max_volume: f64,
    /// Smallest accepted fraction of indexed peaks.
    pub min_completeness: f64,
    /// Largest number of solved assignments tried per family.
    pub max_evaluations: usize,
}

impl Default for IndexingOptions {
    fn default() -> Self {
        Self {
            systems: vec![
                LatticeSystem::Cubic,
                LatticeSystem::Tetragonal,
                LatticeSystem::Hexagonal,
                LatticeSystem::Orthorhombic,
                LatticeSystem::Monoclinic,
            ],
            max_index: 3,
            tolerance: 0.05,
            max_volume: 20_000.0,
            min_completeness: 0.9,
            max_evaluations: 20_000,
        }
    }
}

/// Ranked indexing result.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IndexingResult {
    /// Candidates, best first.
    pub candidates: Vec<IndexedCell>,
}

impl IndexingResult {
    /// Best candidate, when one passed the completeness threshold.
    #[must_use]
    pub fn best(&self) -> Option<&IndexedCell> {
        self.candidates.first()
    }
}

/// Canonical key of the lattice, independent of the cell setting.
type ReducedKey = [i64; 6];

/// One step of the key is 0.01 Angstrom or 0.1 degree; a candidate within one
/// step of a kept key describes the same lattice and is merged with it.
fn keys_close(first: &ReducedKey, second: &ReducedKey) -> bool {
    first
        .iter()
        .zip(second.iter())
        .all(|(a, b)| (a - b).abs() <= 1)
}

fn reduced_key(cell: &CellMetric) -> ReducedKey {
    let Ok(lattice) =
        Lattice::from_parameters(cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)
    else {
        return cell.key();
    };
    let Ok((reduced, _)) = reduce_lattice(&lattice, Reduction::Niggli) else {
        return cell.key();
    };
    let mut parameters = reduced.parameters();
    parameters[..3].sort_by(|first, second| first.total_cmp(second));
    parameters[3..].sort_by(|first, second| first.total_cmp(second));
    [
        (parameters[0] * 100.0).round() as i64,
        (parameters[1] * 100.0).round() as i64,
        (parameters[2] * 100.0).round() as i64,
        (parameters[3] * 10.0).round() as i64,
        (parameters[4] * 10.0).round() as i64,
        (parameters[5] * 10.0).round() as i64,
    ]
}

/// A few unindexed lines are an impurity; 90% counts as fully indexed.
fn completeness_bucket(cell: &IndexedCell) -> i64 {
    (cell.completeness.min(0.90) * 100.0).round() as i64
}

/// Completeness first, then symmetry, then M20.
///
/// A wrong low-symmetry metric can always add parameters until it indexes the
/// same lines, so among cells with the same coverage the highest symmetry is
/// the physical answer; M20 orders the candidates inside one family.
fn rank(first: &IndexedCell, second: &IndexedCell) -> std::cmp::Ordering {
    completeness_bucket(second)
        .cmp(&completeness_bucket(first))
        .then_with(|| first.system.cmp(&second.system))
        .then_with(|| second.m20.total_cmp(&first.m20))
        .then_with(|| first.cell.volume().total_cmp(&second.cell.volume()))
        .then_with(|| second.f_n.total_cmp(&first.f_n))
}

/// Indexes a peak list over the selected metric families.
///
/// Each family is solved from its first free metric parameters (cubic from one
/// peak, monoclinic from four) and validated against every observed peak. The
/// candidates are ranked by M20, and a cell that several families describe
/// equally well is reported in the highest-symmetry family.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for fewer than two peaks, a non-positive
/// wavelength, an empty family list or invalid options.
pub fn index_pattern(
    peaks: &[Peak],
    wavelength: f64,
    options: &IndexingOptions,
) -> Result<IndexingResult, Error> {
    if peaks.len() < 2 {
        return Err(Error::InvalidInput(
            "indexing needs at least two peaks".to_owned(),
        ));
    }
    if !wavelength.is_finite() || wavelength <= 0.0 {
        return Err(Error::InvalidInput(
            "wavelength must be positive and finite".to_owned(),
        ));
    }
    if options.systems.is_empty() {
        return Err(Error::InvalidInput(
            "indexing needs at least one metric family".to_owned(),
        ));
    }
    if options.max_index < 1
        || !options.tolerance.is_finite()
        || options.tolerance <= 0.0
        || !options.max_volume.is_finite()
        || options.max_volume <= 0.0
        || !(0.0..=1.0).contains(&options.min_completeness)
    {
        return Err(Error::InvalidInput(
            "indexing options are invalid".to_owned(),
        ));
    }
    let mut found = general::index_general(peaks, wavelength, options)?;
    // A supercell can index the same lines with a figure of merit an order of
    // magnitude worse; keeping only the good cells of each coverage class is
    // what lets the symmetry comparison below mean something.
    let mut best: BTreeMap<i64, f64> = BTreeMap::new();
    for cell in &found {
        let bucket = completeness_bucket(cell);
        let entry = best.entry(bucket).or_insert(f64::MIN);
        *entry = entry.max(cell.m20);
    }
    found.retain(|cell| {
        best.get(&completeness_bucket(cell))
            .is_none_or(|value| cell.m20 >= 0.1 * value)
    });
    found.sort_by(rank);
    let mut merged: Vec<(ReducedKey, IndexedCell)> = Vec::new();
    for candidate in found {
        let key = reduced_key(&candidate.cell);
        match merged.iter_mut().find(|(kept, _)| keys_close(kept, &key)) {
            Some((_, kept)) => {
                if candidate.system < kept.system {
                    *kept = candidate;
                }
            }
            None => merged.push((key, candidate)),
        }
    }
    let mut merged = merged
        .into_iter()
        .map(|(_, candidate)| candidate)
        .collect::<Vec<_>>();
    merged.sort_by(rank);
    Ok(IndexingResult { candidates: merged })
}
