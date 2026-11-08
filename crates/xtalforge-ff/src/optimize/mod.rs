//! Geometry optimization over the asymmetric unit (Tasks 4.4.1 to 4.4.5).
//!
//! [`fire`] and [`lbfgs`] minimize any [`Calculator`] over the variables
//! [`dof`](self) describes: the sites of the asymmetric unit, regenerated into
//! a P1 cell at every evaluation, plus six symmetric strain components when
//! `cell` is on.  Forces are carried back through the symmetry operations and
//! projected onto each site's allowed displacement, so a special position
//! stays special and a frozen site stays frozen.
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
mod lbfgs;
mod tensor;

pub mod eos;

use serde::{Deserialize, Serialize};
use xtalforge_core::Structure;

use crate::Error;
use crate::calculator::Calculator;
use crate::structure_util::with_frac;

pub use eos::{EosFit, EosForm, MINIMUM_EOS_POINTS, birch_murnaghan, fit_eos, murnaghan};
pub use fire::fire;
pub use lbfgs::lbfgs;

use dof::Dof;

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
        }
    }
}

impl OptimizeOptions {
    /// Checks that the options are usable.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a zero step limit, a non-positive
    /// force or step tolerance, or a negative energy tolerance.
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
        Ok(())
    }
}

/// What an optimization produced.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct OptimizationResult {
    /// Optimized structure, in the same representation as the input.
    pub structure: Structure,
    /// Final energy in eV.
    pub energy: f64,
    /// Largest force on a free site in eV/Angstrom.
    pub max_force: f64,
    /// Number of outer iterations performed.
    pub steps: usize,
    /// Whether the force and energy criteria were met.
    pub converged: bool,
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
    match minimizer {
        Minimizer::Fire => fire(calculator, structure, options),
        Minimizer::Lbfgs => lbfgs(calculator, structure, options),
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

/// One energy and gradient evaluation in the optimizer's variables.
#[derive(Debug, Clone)]
pub(crate) struct Evaluation {
    pub energy: f64,
    pub gradient: Vec<f64>,
    pub max_force: f64,
    pub cell_force: f64,
    pub structure: Structure,
}

/// The energy landscape the minimizers walk.
pub(crate) struct Problem<'a> {
    calculator: &'a dyn Calculator,
    dof: Dof,
}

impl<'a> Problem<'a> {
    pub(crate) fn new(
        calculator: &'a dyn Calculator,
        structure: &Structure,
        options: &OptimizeOptions,
    ) -> Result<Self, Error> {
        options.validate()?;
        let dof = Dof::new(structure, options.mode, options.cell, &options.fixed_atoms)?;
        Ok(Self { calculator, dof })
    }

    pub(crate) fn start(&self) -> Vec<f64> {
        self.dof.start()
    }

    pub(crate) fn project(&self, direction: &mut [f64]) {
        self.dof.project_direction(direction);
    }

    pub(crate) fn evaluate(&mut self, variables: &[f64]) -> Result<Evaluation, Error> {
        let structure = self.dof.structure_of(variables)?;
        let (energy, forces) = self.calculator.energy_forces(&structure)?;
        let site_gradient = self.dof.site_gradient(&forces)?;
        let max_force = self.dof.max_force(&site_gradient);
        let strain = if self.dof.cell() {
            Some(self.strain_gradient(variables, &structure)?)
        } else {
            None
        };
        let gradient = self.dof.flatten(&site_gradient, strain);
        let cell_force = cell_force(&gradient, self.dof.n_variables(), self.dof.cell());
        Ok(Evaluation {
            energy,
            gradient,
            max_force,
            cell_force,
            structure,
        })
    }

    /// Strain gradient from the calculator's stress, or by central differences
    /// when the engine does not provide one.
    fn strain_gradient(
        &mut self,
        variables: &[f64],
        structure: &Structure,
    ) -> Result<[f64; 6], Error> {
        if let Some(stress) = self.calculator.stress(structure)? {
            let mut converted = [[0.0; 3]; 3];
            for (row, values) in stress.iter().enumerate() {
                for (column, value) in values.iter().enumerate() {
                    if let Some(slot) = converted
                        .get_mut(row)
                        .and_then(|values| values.get_mut(column))
                    {
                        *slot = value * GPA_TO_EV_PER_A3;
                    }
                }
            }
            return Ok(self.dof.strain_gradient(variables, converted));
        }
        let base = self.dof.n_variables() - 6;
        let step = 1e-4;
        let mut gradient = [0.0; 6];
        for (index, slot) in gradient.iter_mut().enumerate() {
            let mut high = variables.to_vec();
            let mut low = variables.to_vec();
            if let Some(value) = high.get_mut(base + index) {
                *value += step;
            }
            if let Some(value) = low.get_mut(base + index) {
                *value -= step;
            }
            let energy_high = self.calculator.energy(&self.dof.structure_of(&high)?)?;
            let energy_low = self.calculator.energy(&self.dof.structure_of(&low)?)?;
            *slot = (energy_high - energy_low) / (2.0 * step);
        }
        Ok(gradient)
    }

    pub(crate) fn converged(
        &self,
        evaluation: &Evaluation,
        previous_energy: Option<f64>,
        options: &OptimizeOptions,
    ) -> bool {
        if evaluation.max_force > options.force_tolerance {
            return false;
        }
        if !self.dof.cell() {
            return true;
        }
        if evaluation.cell_force <= options.force_tolerance {
            return true;
        }
        previous_energy
            .map(|energy| (evaluation.energy - energy).abs() <= options.energy_tolerance)
            .unwrap_or(false)
    }

    pub(crate) fn result(
        &self,
        evaluation: Evaluation,
        steps: usize,
        converged: bool,
    ) -> OptimizationResult {
        OptimizationResult {
            structure: evaluation.structure,
            energy: evaluation.energy,
            max_force: evaluation.max_force,
            steps,
            converged,
        }
    }

    pub(crate) fn cap_step(&self, direction: &mut [f64], max_step: f64) {
        let mut longest: f64 = 0.0;
        let sites = self
            .dof
            .n_variables()
            .saturating_sub(if self.dof.cell() { 6 } else { 0 });
        let mut index = 0;
        while index + 2 < sites {
            let vector = [
                direction.get(index).copied().unwrap_or(0.0),
                direction.get(index + 1).copied().unwrap_or(0.0),
                direction.get(index + 2).copied().unwrap_or(0.0),
            ];
            longest = longest.max(crate::math::norm(vector));
            index += 3;
        }
        if longest > max_step {
            let factor = max_step / longest;
            for value in direction.iter_mut() {
                *value *= factor;
            }
        }
    }
}

/// Largest strain gradient component per atom, in eV/Angstrom.
fn cell_force(gradient: &[f64], n_variables: usize, cell: bool) -> f64 {
    if !cell || n_variables < 6 {
        return 0.0;
    }
    let base = n_variables - 6;
    let mut maximum: f64 = 0.0;
    for index in base..n_variables {
        maximum = maximum.max(gradient.get(index).copied().unwrap_or(0.0).abs());
    }
    maximum
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
