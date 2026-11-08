//! Whole-range and structural refinement (Tasks 7.3.2, 7.8 and 7.9.1).
//!
//! [`refine_phase`] builds the calculated pattern of one phase, fits the scale
//! and a constant, linear or quadratic background **analytically**, then runs a
//! small bounded Gauss-Newton loop on the requested structural parameters:
//!
//! * fractional coordinates of every asymmetric site (three per site),
//! * one global isotropic ADP multiplier applied to every site's `U_iso`,
//! * one occupancy per asymmetric site.
//!
//! The structure is edited through [`StructureBuilder`] with the original
//! Hall symbol, lattice and symmetry operations, so the space group is kept and
//! the result is the input structure with new numbers on the same sites. Atoms
//! are never added, removed or bonded.
//!
//! # Approximations
//!
//! * The least-squares step uses central finite differences of the *current*
//!   calculated pattern; the scale and background are refitted analytically
//!   after every trial step, but their derivatives are not part of the
//!   Jacobian.
//! * Position steps are limited to 0.05 fractional units per iteration, ADP
//!   scale steps to 0.25 and occupancy steps to 0.10, so a single iteration can
//!   never move an atom across the cell.
//! * Sites without a `U_iso` start from a 0.005 Angstrom^2 baseline that the
//!   global multiplier scales; sites on special positions should not free
//!   coordinates, because the builder keeps the operations rather than deriving
//!   Wyckoff-allowed directions (Task 2.4.3 does that separately).
//! * The scattering factor, profile and Lorentz-polarisation approximations are
//!   those documented on [`crate::pxrd`]; Rwp inherits them. The model uses the
//!   same area-normalised Gaussian as the PXRD engine, so
//!   [`RietveldOptions::kind`] does not change the profile here.
//! * [`RefineOptions::combined_energy`] is the weight `w` of the combined
//!   objective. This signature has no force-field calculator, so the energy term
//!   is a constant zero during the fit; the weight is used for the reported
//!   `combined_objective` entry, and callers with an energy model should pass
//!   its normalised value to [`combined_objective`] themselves.

mod fit;
mod model;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};
use xtalforge_io::xy::PatternPoint;

use self::fit::{gradient, normal_matrix, solve_linear, weighted_residual};
use self::model::{ModelContext, bounds, build_parameters, clamp, initial_value, max_step};
use crate::peaks::{PeakOptions, find_peaks};
use crate::profiles::Background;
use crate::rietveld::RietveldOptions;

/// Options for [`refine_phase`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefineOptions {
    /// Maximum number of outer least-squares iterations.
    pub max_iterations: usize,
    /// Relative convergence tolerance on the Rwp improvement and step size.
    pub tolerance: f64,
    /// Refine the fractional coordinates of every asymmetric site.
    pub refine_positions: bool,
    /// Refine a global isotropic ADP multiplier.
    pub refine_adps: bool,
    /// Refine one occupancy per asymmetric site.
    pub refine_occupancies: bool,
    /// Optional weight `w` for the combined energy-plus-fit objective.
    pub combined_energy: Option<f64>,
}

impl Default for RefineOptions {
    fn default() -> Self {
        Self {
            max_iterations: 20,
            tolerance: 1e-4,
            refine_positions: false,
            refine_adps: false,
            refine_occupancies: false,
            combined_energy: None,
        }
    }
}

/// A refined phase: the input structure with new numbers on the same sites.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RefinedPhase {
    /// Structure after refinement.
    pub structure: Structure,
    /// Weighted profile R factor of the final model.
    pub rwp: f64,
    /// Refined parameters in a deterministic order: scale, background terms,
    /// then the requested structural parameters.
    pub parameters: Vec<(String, f64)>,
    /// Whether the fit reached a stationary point within the tolerance.
    pub converged: bool,
    /// Outer iterations used.
    pub iterations: usize,
}

/// Combined energy-plus-fit objective of Task 7.9.1.
///
/// The documented objective is
///
/// ```text
/// f = (1 - w) * chi^2 / chi^2_0 + w * (E - E_0) / Delta E
/// ```
///
/// This function expects both arguments **already normalised** by their
/// starting values (`chi^2 / chi^2_0` and `(E - E_0) / Delta E`), so a weight
/// means the same thing on a small cell and on a large framework. The weight is
/// clamped to `[0, 1]`. A non-finite argument returns `f64::INFINITY` rather
/// than panicking.
#[must_use]
pub fn combined_objective(energy: f64, rwp: f64, weight: f64) -> f64 {
    if !energy.is_finite() || !rwp.is_finite() {
        return f64::INFINITY;
    }
    let weight = if weight.is_finite() {
        weight.clamp(0.0, 1.0)
    } else {
        0.0
    };
    (1.0 - weight) * rwp + weight * energy
}

/// Refines one phase against an observed pattern.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for fewer than three points, non-finite or
/// non-monotone observed data, invalid options, a pattern without intensity,
/// a pattern with fewer peaks than the requested structural parameters, or a
/// degenerate scale/background fit.
pub fn refine_phase(
    structure: &Structure,
    observed: &[PatternPoint],
    options: &RietveldOptions,
    refine: &RefineOptions,
) -> Result<RefinedPhase, Error> {
    validate_observed(observed)?;
    if !refine.tolerance.is_finite() || refine.tolerance <= 0.0 {
        return Err(Error::InvalidInput(
            "refinement tolerance must be positive and finite".to_owned(),
        ));
    }
    if let Some(weight) = refine.combined_energy
        && (!weight.is_finite() || !(0.0..=1.0).contains(&weight))
    {
        return Err(Error::InvalidInput(format!(
            "combined energy weight must be in [0, 1], got {weight}"
        )));
    }
    let terms = background_terms(options.background);
    let parameters = build_parameters(structure, refine);
    let peaks = count_peaks(observed)?;
    if parameters.len() > peaks {
        return Err(Error::InvalidInput(format!(
            "refinement would free {} structural parameters but the observed pattern supports only {peaks} peaks; disable positions, ADPs or occupancies, or observe a wider range",
            parameters.len()
        )));
    }
    let base_sites = structure.sites().to_vec();
    let base_u: Vec<f64> = base_sites
        .iter()
        .map(|site| site.u_iso().unwrap_or(model::DEFAULT_U_ISO))
        .collect();
    let mut values: Vec<f64> = parameters
        .iter()
        .map(|parameter| initial_value(parameter.kind, &base_sites))
        .collect();
    let pxrd = model::pxrd_options(options, observed);
    let context = ModelContext {
        base: structure,
        base_u: &base_u,
        parameters: &parameters,
        observed,
        pxrd: &pxrd,
    };
    let mut current = context.model(&values, terms)?;
    let initial_rwp = current.rwp.max(1e-12);
    let mut iterations = 0;
    let mut converged = parameters.is_empty();
    let mut lambda = 1e-3;

    while !parameters.is_empty() && iterations < refine.max_iterations && !converged {
        iterations += 1;
        let residual = weighted_residual(observed, &current.model);
        let jacobian = context.jacobian(&values)?;
        let normal = normal_matrix(&jacobian);
        let gradient = gradient(&jacobian, &residual);
        let rhs: Vec<f64> = gradient.iter().map(|value| -value).collect();
        let mut accepted = false;
        for _ in 0..6 {
            let mut damped = normal.clone();
            for index in 0..parameters.len() {
                if let Some(slot) = damped.get_mut(index).and_then(|row| row.get_mut(index)) {
                    *slot += lambda * (1.0 + slot.abs());
                }
            }
            let Some(step) = solve_linear(&damped, &rhs) else {
                lambda = (lambda * 10.0).min(1e12);
                continue;
            };
            let mut candidate_values = values.clone();
            let mut largest = 0.0_f64;
            for (index, parameter) in parameters.iter().enumerate() {
                let (low, high) = bounds(parameter.kind);
                let limit = max_step(parameter.kind);
                let raw = step.get(index).copied().unwrap_or(0.0);
                let delta = raw.clamp(-limit, limit);
                let Some(value) = values.get(index).copied() else {
                    continue;
                };
                let candidate = clamp(value + delta, low, high);
                if let Some(slot) = candidate_values.get_mut(index) {
                    *slot = candidate;
                }
                largest = largest.max(delta.abs() / (1.0 + value.abs()));
            }
            if largest < refine.tolerance {
                converged = true;
                break;
            }
            let candidate = context.model(&candidate_values, terms)?;
            if candidate.rwp < current.rwp {
                let improvement = (current.rwp - candidate.rwp) / initial_rwp;
                values = candidate_values;
                current = candidate;
                lambda = (lambda / 3.0).max(1e-12);
                accepted = true;
                if improvement < refine.tolerance {
                    converged = true;
                }
                break;
            }
            lambda = (lambda * 5.0).min(1e12);
        }
        if !accepted {
            converged = true;
        }
    }

    let mut reported = vec![("scale".to_owned(), current.scale)];
    for (index, coefficient) in current.background.iter().enumerate() {
        reported.push((format!("background_{index}"), *coefficient));
    }
    for (parameter, value) in parameters.iter().zip(values.iter()) {
        reported.push((parameter.name.clone(), *value));
    }
    if let Some(weight) = refine.combined_energy {
        let merit = combined_objective(0.0, current.rwp / initial_rwp, weight);
        reported.push(("combined_objective".to_owned(), merit));
    }
    Ok(RefinedPhase {
        structure: current.structure,
        rwp: current.rwp,
        parameters: reported,
        converged,
        iterations,
    })
}

fn validate_observed(observed: &[PatternPoint]) -> Result<(), Error> {
    if observed.len() < 3 {
        return Err(Error::InvalidInput(
            "refinement needs at least three observed points".to_owned(),
        ));
    }
    for (index, point) in observed.iter().enumerate() {
        if !point.two_theta.is_finite() || !point.intensity.is_finite() {
            return Err(Error::InvalidInput(format!(
                "observed point {index} must be finite"
            )));
        }
        if let Some(previous) = index.checked_sub(1).and_then(|value| observed.get(value))
            && point.two_theta <= previous.two_theta
        {
            return Err(Error::InvalidInput(format!(
                "observed 2 theta must be strictly increasing at point {index}"
            )));
        }
    }
    Ok(())
}

fn background_terms(background: Background) -> usize {
    match background {
        Background::Constant(_) => 1,
        Background::Linear(..) => 2,
        Background::Quadratic(..) => 3,
    }
}

fn count_peaks(observed: &[PatternPoint]) -> Result<usize, Error> {
    let mut steps = Vec::new();
    for pair in observed.windows(2) {
        if let [first, second] = pair
            && second.two_theta > first.two_theta
        {
            steps.push(second.two_theta - first.two_theta);
        }
    }
    steps.sort_by(f64::total_cmp);
    let step = steps
        .get(steps.len() / 2)
        .copied()
        .unwrap_or(0.02)
        .max(0.001);
    let peaks = find_peaks(
        observed,
        &PeakOptions {
            min_height: 0.0,
            min_distance: (2.0 * step).max(0.01),
            smoothing: 1,
        },
    )?;
    Ok(peaks.len())
}
