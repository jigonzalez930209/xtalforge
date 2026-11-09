//! Pawley whole-pattern fit with esds and free cell/instrument parameters
//! (Task 7.6.2).
//!
//! The observed pattern is modelled as
//!
//! `y(x) = sum_k I_k P(x; 2theta_k, w_k) + sum_j c_j T_j(x)`
//!
//! where `P` is the profile, `I_k` the intensity of reflection `k`, `w_k` the
//! width from the size/strain model and `T_j` a Chebyshev background. For a
//! fixed cell, zero shift and width model the model is **linear in the
//! intensities and background coefficients**, so one weighted least-squares
//! solve gives them together with their covariance.
//!
//! The cell, zero shift and width model are non-linear because they move the
//! reflection positions. They are refined by alternating: solve the linear
//! part, refine the non-linear parameters with those intensities fixed, and
//! solve again. This is the standard Pawley separation and keeps the outer
//! problem small.
//!
//! **Documented conventions.** Intensities are profile heights, not
//! integrated areas, and they are reported per powder line (multiplicity
//! included). The size term is a constant FWHM in degrees; the strain term
//! adds `strain * tan(theta)` in quadrature. Weights are `1 / uncertainty^2`
//! when the point carries one, otherwise `1 / max(y, 1)`. esds come from the
//! diagonal of the covariance at the solution, scaled by the reduced
//! chi-squared; the non-linear esds use a finite-difference Jacobian with the
//! intensities fixed at their fitted values and are therefore approximate.
//!
//! A Bragg-Brentano cell also carries a systematic error of order `1e-4`
//! Angstrom from the sample displacement that no esd includes; the reference
//! says the same and this note is the documented deviation.

#![allow(clippy::indexing_slicing)] // bounded local matrices

mod cell;
mod esd;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};
use xtalforge_io::xy::PatternPoint;

use crate::profiles::ProfileKind;
use crate::pxrd::{PatternPoint as PatternValue, PxrdOptions, Radiation, calculate_pxrd};

use self::cell::Family;

/// Options for the Pawley fit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PawleyOptions {
    /// Profile shape.
    pub kind: ProfileKind,
    /// Pseudo-Voigt mixing.
    pub eta: f64,
    /// Initial FWHM in degrees (also the fixed width when size/strain are off).
    pub fwhm: f64,
    /// Chebyshev background order (0 is a constant).
    pub background_order: usize,
    /// Refine the free cell parameters.
    pub refine_cell: bool,
    /// Refine a zero shift in degrees.
    pub refine_zero: bool,
    /// Refine the constant FWHM (size) term.
    pub refine_size: bool,
    /// Refine the `tan(theta)` (strain) term.
    pub refine_strain: bool,
    /// Iterations of the non-linear step per alternation.
    pub max_iterations: usize,
    /// Largest number of allowed reflections the fit accepts.
    pub max_reflections: usize,
}

impl Default for PawleyOptions {
    fn default() -> Self {
        Self {
            kind: ProfileKind::PseudoVoigt,
            eta: 0.5,
            fwhm: 0.1,
            background_order: 1,
            refine_cell: true,
            refine_zero: true,
            refine_size: false,
            refine_strain: false,
            max_iterations: 20,
            max_reflections: 400,
        }
    }
}

impl PawleyOptions {
    fn validate(&self) -> Result<(), Error> {
        if !self.fwhm.is_finite() || self.fwhm <= 0.0 {
            return Err(Error::InvalidInput(
                "the Pawley FWHM must be positive and finite".to_owned(),
            ));
        }
        if !self.eta.is_finite() || !(0.0..=1.0).contains(&self.eta) {
            return Err(Error::InvalidInput(
                "the Pawley mixing must be in [0, 1]".to_owned(),
            ));
        }
        if self.max_reflections == 0 {
            return Err(Error::InvalidInput(
                "the Pawley fit needs room for at least one reflection".to_owned(),
            ));
        }
        if self.background_order > 12 {
            return Err(Error::InvalidInput(
                "the Pawley background order must not exceed 12".to_owned(),
            ));
        }
        Ok(())
    }
}

/// One fitted reflection.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PawleyReflection {
    /// Miller indices.
    pub hkl: [i32; 3],
    /// Interplanar spacing in Angstrom.
    pub d: f64,
    /// Bragg position in degrees two-theta.
    pub two_theta: f64,
    /// Fitted intensity (profile height).
    pub intensity: f64,
    /// Standard uncertainty of `intensity`.
    pub esd: f64,
}

/// Result of a Pawley fit.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PawleyResult {
    /// Fitted reflections, ordered by two-theta.
    pub reflections: Vec<PawleyReflection>,
    /// Cell parameters `(a, b, c, alpha, beta, gamma)` after refinement.
    pub cell: [f64; 6],
    /// Standard uncertainties of the cell parameters (zero when fixed).
    pub cell_esds: [f64; 6],
    /// Zero shift in degrees.
    pub zero_shift: f64,
    /// Standard uncertainty of the zero shift.
    pub zero_shift_esd: f64,
    /// Constant FWHM (size) term in degrees.
    pub size: f64,
    /// `tan(theta)` (strain) term in degrees.
    pub strain: f64,
    /// Profile R factor.
    pub rp: f64,
    /// Weighted profile R factor.
    pub rwp: f64,
    /// Expected R factor.
    pub rexp: f64,
    /// Reduced chi-squared.
    pub chi_squared: f64,
    /// Calculated pattern at the observed points.
    pub pattern: Vec<PatternValue>,
    /// Whether the non-linear step reported convergence.
    pub converged: bool,
    /// Alternations performed.
    pub iterations: usize,
}

/// Fits a Pawley model to a measured pattern.
///
/// The reflections allowed by the structure's space group are taken from the
/// PXRD engine; their positions follow the refined cell and zero shift, and
/// their intensities are free parameters.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for fewer than five points, non-monotone
/// data, invalid options, no reflection in range or a singular linear system.
pub fn pawley_refine(
    structure: &Structure,
    observed: &[PatternPoint],
    options: &PawleyOptions,
) -> Result<PawleyResult, Error> {
    options.validate()?;
    validate_observed(observed)?;
    let base = structure.lattice().parameters();
    let family = Family::from_parameters(base);
    let first = observed.first().map_or(0.0, |point| point.two_theta);
    let last = observed.last().map_or(0.0, |point| point.two_theta);
    let pattern = calculate_pxrd(
        structure,
        &PxrdOptions {
            radiation: Radiation::CuKa,
            two_theta_range: [first.max(0.0), last.min(180.0)],
            step: 0.02,
            fwhm: options.fwhm,
        },
    )?;
    let mut hkls: Vec<[i32; 3]> = Vec::new();
    for reflection in &pattern.reflections {
        if reflection.intensity > 1e-8 && !hkls.contains(&reflection.hkl) {
            hkls.push(reflection.hkl);
        }
    }
    if hkls.is_empty() {
        return Err(Error::InvalidInput(
            "the structure has no allowed reflection in the observed range".to_owned(),
        ));
    }
    if hkls.len() > options.max_reflections {
        return Err(Error::InvalidInput(format!(
            "the cell predicts {} reflections; the Pawley fit is limited to {}",
            hkls.len(),
            options.max_reflections
        )));
    }
    let wavelength = Radiation::CuKa.wavelength()?;
    let layout = Layout::new(family, base, options);
    let mut parameters = layout.initial(options);
    let mut converged = false;
    let mut iterations = 0;
    let mut design = Design::new(observed, &hkls, options, &layout, &parameters, wavelength)?;
    let mut linear = design.solve()?;
    if !parameters.is_empty() && options.max_iterations > 0 {
        let steps = layout.steps();
        let (refined, _, reached) =
            nelder_mead(&parameters, &steps, |values| design.chi_squared(values))?;
        parameters = refined;
        converged = reached;
        iterations = 1;
        design = design.refit(&parameters)?;
        linear = design.solve()?;
    }
    let parameter_count = design.parameter_count() + layout.len();
    let variance = design.reduced_variance(&parameters);
    let intensities = (0..linear.intensities.len())
        .map(|index| (linear.value(index), linear.esd(index, variance)))
        .collect::<Vec<_>>();
    let (cell_esds, zero_esd) = layout.errors(&parameters, &design, &linear, variance);
    let fit_values = design.model(&linear);
    let (rp, rwp, rexp, chi_squared) = factors(observed, &fit_values, parameter_count);
    let reflections = design
        .reflections
        .iter()
        .zip(intensities.iter())
        .map(|(reflection, (intensity, esd))| PawleyReflection {
            hkl: reflection.hkl,
            d: reflection.d,
            two_theta: reflection.two_theta,
            intensity: *intensity,
            esd: *esd,
        })
        .collect();
    Ok(PawleyResult {
        reflections,
        cell: layout.cell(&parameters),
        cell_esds,
        zero_shift: layout.zero(&parameters),
        zero_shift_esd: zero_esd,
        size: layout.size_value(&parameters),
        strain: layout.strain_value(&parameters),
        rp,
        rwp,
        rexp,
        chi_squared,
        pattern: observed
            .iter()
            .zip(fit_values.iter())
            .map(|(point, value)| PatternValue {
                two_theta: point.two_theta,
                intensity: *value,
            })
            .collect(),
        converged,
        iterations,
    })
}

fn validate_observed(observed: &[PatternPoint]) -> Result<(), Error> {
    if observed.len() < 5 {
        return Err(Error::InvalidInput(
            "the Pawley fit needs at least five observed points".to_owned(),
        ));
    }
    for (index, point) in observed.iter().enumerate() {
        if !point.two_theta.is_finite() || !point.intensity.is_finite() {
            return Err(Error::InvalidInput(format!(
                "observed point {index} must be finite"
            )));
        }
        if let Some(previous) = index.checked_sub(1).and_then(|slot| observed.get(slot))
            && point.two_theta <= previous.two_theta
        {
            return Err(Error::InvalidInput(format!(
                "observed 2 theta must be strictly increasing at point {index}"
            )));
        }
    }
    Ok(())
}

mod design;
mod model;
mod optimize;

use self::design::Design;
use self::model::{Layout, factors};
use self::optimize::nelder_mead;
