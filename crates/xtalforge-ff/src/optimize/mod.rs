//! Geometry optimization over the asymmetric unit (Tasks 4.4.1 to 4.4.5,
//! extended by Task 14.2).
//!
//! [`fire`] and [`lbfgs`] minimize any [`Calculator`] over the variables
//! [`dof`](self) describes: the sites of the asymmetric unit, regenerated into
//! a P1 cell at every evaluation, plus six symmetric strain components when
//! `cell` is on.  Forces are carried back through the symmetry operations and
//! projected onto each site's allowed displacement, so a special position
//! stays special and a frozen site stays frozen.
//!
//! [`CellFreedom`] narrows the cell variables: held lengths, angles and volume
//! are re-imposed after every cell update (Task 14.2.1).  `pressure_gpa` adds
//! the `P V` term to the objective, so a variable-cell run at pressure
//! minimizes the enthalpy while the result still reports the potential energy.
//!
//! [`OptimizeOptions::trajectory`] records frames, [`RunControl`] reports
//! progress and stops the run, and the result carries both the frames and a
//! `cancelled` flag (Task 14.2.2).  The older [`fire`], [`lbfgs`] and
//! [`minimize`] entry points keep their signature and behave as before.
//!
//! The result contains a ready structure.  [`OptimizationResult::apply_to`]
//! writes the optimized coordinates and cell back onto the input structure,
//! preserving its labels, charges, bonds and metadata (Task 4.4.4).
//!
//! Energies at this boundary are eV, forces eV/Angstrom and the maximum force
//! is over the free sites.  The cell variables are compared with the same
//! tolerance after dividing their gradient by the number of atoms, so one
//! number means one thing whether a cell is relaxing or not.

mod dof;
mod fire;
mod freedom;
mod lbfgs;
mod problem;
mod symmetry;
mod tensor;

pub mod eos;

use serde::{Deserialize, Serialize};
use xtalforge_core::Structure;

use crate::Error;
use crate::calculator::Calculator;
use crate::run::RunControl;
use crate::structure_util::with_frac;
use crate::trajectory::{Trajectory, TrajectoryOptions};

pub use eos::{EosFit, EosForm, MINIMUM_EOS_POINTS, birch_murnaghan, fit_eos, murnaghan};
pub use fire::{fire, fire_with};
pub use freedom::CellFreedom;
pub use lbfgs::{lbfgs, lbfgs_with};

use problem::{Evaluation, Problem};

/// `1 GPa` in `eV/Angstrom^3`.
pub(crate) const GPA_TO_EV_PER_A3: f64 = 1.0 / 160.217_663_4;

/// Which variables the optimizer is allowed to move.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum RelaxMode {
    /// Every site is free; the structure must be P1.
    P1,
    /// The asymmetric unit is the variable and the space group is preserved.
    #[default]
    Asymmetric,
}

/// Options for one geometry optimization.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OptimizeOptions {
    /// Maximum number of outer iterations.
    pub max_steps: usize,
    /// Largest allowed force on a free site, in eV/Angstrom.
    pub force_tolerance: f64,
    /// Energy change below which a run may call itself converged, in eV.
    pub energy_tolerance: f64,
    /// Whether the six symmetric strain components are variables.
    pub cell: bool,
    /// Indices of the asymmetric sites that must not move.
    pub fixed_atoms: Vec<usize>,
    /// Variables to optimize.
    pub mode: RelaxMode,
    /// Largest Cartesian displacement of one site in one step, in Angstrom.
    pub max_step: f64,
    /// Which cell lengths, angles and volume may change when `cell` is on.
    ///
    /// The default, [`CellFreedom::full`], reproduces the behaviour of earlier
    /// releases exactly.
    #[serde(default)]
    pub cell_freedom: CellFreedom,
    /// External pressure in GPa, applied as the `P V` term to the objective.
    ///
    /// Zero (the default) minimizes the potential energy.  The value is only
    /// meaningful with `cell` on, because a fixed cell has a constant volume.
    #[serde(default)]
    pub pressure_gpa: f64,
    /// Trajectory recording cadence; off by default.
    #[serde(default)]
    pub trajectory: TrajectoryOptions,
}

impl Default for OptimizeOptions {
    fn default() -> Self {
        Self {
            max_steps: 500,
            force_tolerance: 0.05,
            energy_tolerance: 1e-6,
            cell: false,
            fixed_atoms: Vec::new(),
            mode: RelaxMode::default(),
            max_step: 0.2,
            cell_freedom: CellFreedom::default(),
            pressure_gpa: 0.0,
            trajectory: TrajectoryOptions::default(),
        }
    }
}

impl OptimizeOptions {
    /// Checks that the options are usable.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a zero step limit, a non-positive
    /// force or step tolerance, a negative energy tolerance, a non-finite
    /// pressure, or a cell freedom that cannot drive a variable-cell run.
    pub fn validate(&self) -> Result<(), Error> {
        if self.max_steps == 0 {
            return Err(Error::InvalidInput(
                "an optimization needs at least one step".to_owned(),
            ));
        }
        if !self.force_tolerance.is_finite() || self.force_tolerance <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "force tolerance must be positive and finite, got {}",
                self.force_tolerance
            )));
        }
        if !self.energy_tolerance.is_finite() || self.energy_tolerance < 0.0 {
            return Err(Error::InvalidInput(format!(
                "energy tolerance must be non-negative and finite, got {}",
                self.energy_tolerance
            )));
        }
        if !self.max_step.is_finite() || self.max_step <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "maximum step must be positive and finite, got {}",
                self.max_step
            )));
        }
        if !self.pressure_gpa.is_finite() {
            return Err(Error::InvalidInput(format!(
                "external pressure must be finite, got {}",
                self.pressure_gpa
            )));
        }
        if self.cell {
            self.cell_freedom.validate()?;
        }
        Ok(())
    }
}

/// What an optimization produced.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OptimizationResult {
    /// Optimized structure, in the same representation as the input.
    pub structure: Structure,
    /// Final potential energy in eV.
    pub energy: f64,
    /// Final optimization objective in eV, `energy + pressure * volume`.
    ///
    /// Equal to `energy` at the default zero pressure.
    pub enthalpy: f64,
    /// Largest force on a free site in eV/Angstrom.
    pub max_force: f64,
    /// Number of outer iterations performed.
    pub steps: usize,
    /// Whether the force and energy criteria were met.
    pub converged: bool,
    /// Whether the run stopped because a [`RunControl`] asked it to.
    pub cancelled: bool,
    /// Recorded frames, empty unless the options enabled a trajectory.
    pub trajectory: Trajectory,
}

impl OptimizationResult {
    /// Writes the optimized coordinates and cell onto `structure`.
    ///
    /// The input keeps its labels, elements, charges, occupancies, bonds and
    /// metadata; only `frac` and the lattice come from this result.  The input
    /// must have the same number of asymmetric sites, which it does when the
    /// optimization was started from it (Task 4.4.4).
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the site counts differ, or a core
    /// error when the structure cannot be rebuilt.
    pub fn apply_to(&self, structure: &Structure) -> Result<Structure, Error> {
        if structure.sites().len() != self.structure.sites().len() {
            return Err(Error::InvalidInput(format!(
                "cannot apply {} optimized sites to a structure with {} sites",
                self.structure.sites().len(),
                structure.sites().len()
            )));
        }
        let frac: Vec<[f64; 3]> = self
            .structure
            .sites()
            .iter()
            .map(xtalforge_core::Site::frac)
            .collect();
        with_frac(structure, *self.structure.lattice(), &frac)
    }

    /// Fractional coordinates of the optimized sites.
    #[must_use]
    pub fn frac(&self) -> Vec<[f64; 3]> {
        self.structure
            .sites()
            .iter()
            .map(xtalforge_core::Site::frac)
            .collect()
    }
}

/// Which minimizer to run.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Minimizer {
    /// FIRE: molecular dynamics with the velocity steered downhill.
    Fire,
    /// L-BFGS: a limited-memory quasi-Newton line search.
    Lbfgs,
}

/// Runs one minimizer.
///
/// # Errors
///
/// Returns the errors of [`fire`] or [`lbfgs`].
pub fn minimize(
    calculator: &dyn Calculator,
    structure: &Structure,
    options: &OptimizeOptions,
    minimizer: Minimizer,
) -> Result<OptimizationResult, Error> {
    minimize_with(
        calculator,
        structure,
        options,
        minimizer,
        &mut RunControl::new(),
    )
}

/// Runs one minimizer with cancellation and progress reporting.
///
/// # Errors
///
/// Returns the errors of [`fire_with`] or [`lbfgs_with`].
pub fn minimize_with(
    calculator: &dyn Calculator,
    structure: &Structure,
    options: &OptimizeOptions,
    minimizer: Minimizer,
    control: &mut RunControl<'_>,
) -> Result<OptimizationResult, Error> {
    match minimizer {
        Minimizer::Fire => fire_with(calculator, structure, options, control),
        Minimizer::Lbfgs => lbfgs_with(calculator, structure, options, control),
    }
}

/// Writes an optimization result back onto its input structure.
///
/// # Errors
///
/// Returns the errors of [`OptimizationResult::apply_to`].
pub fn apply_optimized_geometry(
    result: &OptimizationResult,
    structure: &Structure,
) -> Result<Structure, Error> {
    result.apply_to(structure)
}

/// Records one frame when the cadence asks for it.
pub(crate) fn record(
    trajectory: &mut Trajectory,
    cadence: TrajectoryOptions,
    step: usize,
    point: Option<usize>,
    evaluation: &Evaluation,
) {
    if cadence.enabled() && step.is_multiple_of(cadence.every) {
        trajectory.push(evaluation.frame(step, point, cadence));
    }
}

/// Appends an optimizer trajectory to a scan trajectory.
pub(crate) fn append_frames(target: &mut Trajectory, source: Trajectory, point: usize) {
    for mut frame in source.frames {
        frame.point = Some(point);
        target.push(frame);
    }
}

/// Dot product of two flat vectors.
pub(crate) fn dot(first: &[f64], second: &[f64]) -> f64 {
    first.iter().zip(second).map(|(a, b)| a * b).sum()
}

/// Euclidean norm of a flat vector.
pub(crate) fn norm(vector: &[f64]) -> f64 {
    dot(vector, vector).sqrt()
}

/// `a + factor * b` in place.
pub(crate) fn axpy(target: &mut [f64], factor: f64, source: &[f64]) {
    for (slot, value) in target.iter_mut().zip(source) {
        *slot += factor * value;
    }
}

/// `first - second`.
pub(crate) fn difference(first: &[f64], second: &[f64]) -> Vec<f64> {
    first.iter().zip(second).map(|(a, b)| a - b).collect()
}
