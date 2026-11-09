//! Whole-range and structural refinement (Tasks 7.3.2, 7.8 and 7.9.1).
//!
//! [`refine_phase`] builds the calculated pattern of one phase and refines it
//! in the McCusker order: the scale and background are fitted analytically,
//! then the zero shift and the free cell numbers, then the profile width, and
//! finally the structure (Wyckoff-allowed coordinates, isotropic ADPs,
//! occupancies and March-Dollase texture).
//!
//! The structure is edited through [`xtalforge_core::StructureBuilder`] with the
//! original Hall symbol, lattice and symmetry operations, so the space group is
//! kept and the result is the input structure with new numbers on the same
//! sites. Atoms are never added, removed or bonded.
//!
//! # Approximations
//!
//! * Every stage except the analytic scale/background solve uses a bounded
//!   Nelder-Mead simplex, because the zero shift, cell and width are strongly
//!   correlated with each other and with the scale. The esds come from central
//!   differences of the complete model at the final values, so the analytic
//!   scale/background refit is projected into every Jacobian column.
//! * A site on a special position moves only along the null space of its
//!   stabilizer, computed from the operations; a general position has three
//!   directions, a mirror or axis one, an inversion centre none. The free
//!   subspace is public as [`free_directions`].
//! * Only the numbers the metric family leaves free are refined, and the
//!   family ties (`a = b = c`, `a = b`, right angles) are restored after every
//!   move.
//! * The profile uses the area-normalised Gaussian of the PXRD engine;
//!   `RietveldOptions::kind` is accepted but does not change the curve here.
//!   With [`RefineOptions::doublet`] the copper K alpha 1/K alpha 2 pair is
//!   placed from the K alpha 1 position with the 0.5 intensity ratio.
//! * The ADP is one global multiplier on every site's `U_iso`; a free `U_ij`
//!   per site is not implemented.
//! * [`RefineOptions::combined_energy`] is the weight `w` of the combined
//!   objective. This signature has no force-field calculator, so the energy
//!   term is a constant zero during the fit.

mod curve;
mod fit;
mod model;
mod stage;
mod wyckoff;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure, SymOp};
use xtalforge_io::xy::PatternPoint;

use self::fit::{invert, weighted_residual};
use self::model::{ModelContext, ParameterKind, build_parameters, initial_value, pxrd_options};
use self::stage::run_stage;
use crate::peaks::{PeakOptions, find_peaks};
use crate::profiles::Background;
use crate::rietveld::RietveldOptions;

/// Options for [`refine_phase`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefineOptions {
    /// Maximum number of outer least-squares iterations per stage.
    pub max_iterations: usize,
    /// Relative convergence tolerance on the Rwp improvement and step size.
    pub tolerance: f64,
    /// Refine the Wyckoff-allowed coordinates of every asymmetric site.
    pub refine_positions: bool,
    /// Refine a global isotropic ADP multiplier.
    pub refine_adps: bool,
    /// Refine one occupancy per asymmetric site.
    pub refine_occupancies: bool,
    /// Refine the zero shift.
    pub refine_zero: bool,
    /// Refine the cell numbers the metric family leaves free.
    pub refine_cell: bool,
    /// Refine the profile width multiplier.
    pub refine_width: bool,
    /// Refine the March-Dollase texture coefficient.
    pub refine_texture: bool,
    /// Preferred-orientation axis in fractional direct coordinates.
    pub texture_axis: [f64; 3],
    /// Refine against the copper K alpha 1/K alpha 2 doublet.
    pub doublet: bool,
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
            refine_zero: false,
            refine_cell: false,
            refine_width: false,
            refine_texture: false,
            texture_axis: [0.0, 0.0, 1.0],
            doublet: false,
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
    /// then the requested parameters.
    pub parameters: Vec<(String, f64)>,
    /// Standard uncertainties of the refined parameters.
    pub parameter_esds: Vec<(String, f64)>,
    /// Stages run, in McCusker order.
    pub stages: Vec<String>,
    /// Whether every stage reached a stationary point.
    pub converged: bool,
    /// Outer iterations used by all stages.
    pub iterations: usize,
}

/// Orthonormal directions a site may move along without leaving its Wyckoff
/// position.
///
/// A displacement belongs to the null space of every stabilizer operation:
/// three directions on a general position, one on a mirror or rotation axis,
/// none at an inversion centre. The basis is deterministic and empty when the
/// site is fully fixed.
#[must_use]
pub fn free_directions(operations: &[SymOp], frac: [f64; 3]) -> Vec<[f64; 3]> {
    wyckoff::allowed_directions(operations, frac)
}

/// Combined energy-plus-fit objective of Task 7.9.1.
///
/// This function expects both arguments **already normalised** by their
/// starting values (`chi^2 / chi^2_0` and `(E - E_0) / Delta E`), so a weight
/// means the same thing on a small cell and on a large framework. The weight is
/// clamped to `[0, 1]`. A non-finite argument returns `f64::INFINITY`.
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
    if refine.refine_texture {
        let norm = refine
            .texture_axis
            .iter()
            .map(|value| value * value)
            .sum::<f64>();
        if !norm.is_finite() || norm <= 1e-12 {
            return Err(Error::InvalidInput(
                "the texture axis must be a finite non-zero direction".to_owned(),
            ));
        }
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
            "refinement would free {} parameters but the observed pattern supports only {peaks} peaks; disable positions, ADPs, occupancies, cell or width, or observe a wider range",
            parameters.len()
        )));
    }
    let base_u: Vec<f64> = structure
        .sites()
        .iter()
        .map(|site| site.u_iso().unwrap_or(model::DEFAULT_U_ISO))
        .collect();
    let mut values: Vec<f64> = parameters
        .iter()
        .map(|parameter| initial_value(parameter.kind, structure))
        .collect();
    let pxrd = pxrd_options(options, observed);
    let context = ModelContext {
        base: structure,
        base_u: &base_u,
        parameters: &parameters,
        observed,
        pxrd: &pxrd,
        texture_axis: refine.texture_axis,
        doublet: refine.doublet,
    };
    let mut stages = vec!["background_scale".to_owned()];
    let mut current = context.model(&values, terms)?;
    let mut iterations = 0;
    let mut converged = true;
    type StageSelector = fn(ParameterKind) -> bool;
    let stage_order: [(&str, StageSelector); 3] = [
        ("zero_cell", |kind| {
            matches!(kind, ParameterKind::Zero | ParameterKind::Cell { .. })
        }),
        ("width", |kind| matches!(kind, ParameterKind::WidthScale)),
        ("structure", |kind| {
            matches!(
                kind,
                ParameterKind::Position { .. }
                    | ParameterKind::AdpScale
                    | ParameterKind::Occupancy { .. }
                    | ParameterKind::Texture
            )
        }),
    ];
    for (name, predicate) in stage_order {
        let indices: Vec<usize> = parameters
            .iter()
            .enumerate()
            .filter(|(_, parameter)| predicate(parameter.kind))
            .map(|(index, _)| index)
            .collect();
        if indices.is_empty() {
            continue;
        }
        let result = run_stage(&context, &mut values, terms, &indices, refine)?;
        current = result.model;
        iterations += result.iterations;
        converged &= result.converged;
        stages.push(name.to_owned());
    }

    let mut reported = vec![("scale".to_owned(), current.scale)];
    for (index, coefficient) in current.background.iter().enumerate() {
        reported.push((format!("background_{index}"), *coefficient));
    }
    for (parameter, value) in parameters.iter().zip(values.iter()) {
        reported.push((parameter.name.clone(), *value));
    }
    let parameter_esds = parameter_errors(&context, &values, terms, &parameters)?;
    if let Some(weight) = refine.combined_energy {
        let merit = combined_objective(0.0, current.rwp / current.rwp.max(1e-12), weight);
        reported.push(("combined_objective".to_owned(), merit));
    }
    Ok(RefinedPhase {
        structure: current.structure,
        rwp: current.rwp,
        parameters: reported,
        parameter_esds,
        stages,
        converged,
        iterations,
    })
}

/// Covariance-based esds of the free parameters at the final values.
fn parameter_errors(
    context: &ModelContext,
    values: &[f64],
    terms: usize,
    parameters: &[model::Parameter],
) -> Result<Vec<(String, f64)>, Error> {
    if parameters.is_empty() {
        return Ok(Vec::new());
    }
    let indices: Vec<usize> = (0..parameters.len()).collect();
    let model = context.model(values, terms)?;
    let jacobian = context.jacobian(values, &indices, terms)?;
    let normal = fit::normal_matrix(&jacobian);
    let residual = weighted_residual(context.observed, &model.model);
    let chi_squared: f64 = residual.iter().map(|value| value * value).sum();
    let free = parameters.len() + terms + 1;
    let degrees = context.observed.len().saturating_sub(free).max(1) as f64;
    let variance = chi_squared / degrees;
    let covariance = invert(&normal).unwrap_or_default();
    Ok(parameters
        .iter()
        .enumerate()
        .map(|(index, parameter)| {
            let value = covariance
                .get(index)
                .and_then(|row| row.get(index))
                .copied()
                .unwrap_or(0.0);
            (parameter.name.clone(), (value * variance).max(0.0).sqrt())
        })
        .collect())
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
            ..PeakOptions::default()
        },
    )?;
    Ok(peaks.len())
}
