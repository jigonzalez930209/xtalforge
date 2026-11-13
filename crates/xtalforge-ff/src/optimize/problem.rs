//! The energy landscape the minimizers walk (Task 4.4, extended by Task 14.2).
//!
//! The problem owns the variables ([`Dof`]) and turns a variable vector into an
//! energy, a gradient, the maximum force and the structure the calculator saw.
//! The optimization objective is the enthalpy `E + P V`; at the default zero
//! pressure it is the potential energy, which is what the result reports.

use xtalforge_core::Structure;

use crate::Error;
use crate::calculator::Calculator;
use crate::optimize::dof::Dof;
use crate::optimize::{GPA_TO_EV_PER_A3, OptimizationResult, OptimizeOptions};
use crate::trajectory::{Trajectory, TrajectoryFrame, TrajectoryOptions};

/// One energy and gradient evaluation in the optimizer's variables.
#[derive(Debug, Clone)]
pub(crate) struct Evaluation {
    /// Potential energy in eV.
    pub energy: f64,
    /// Optimization objective in eV, `energy + pressure * volume`.
    pub enthalpy: f64,
    /// Gradient in eV per variable.
    pub gradient: Vec<f64>,
    /// Largest force on a free site in eV/Angstrom.
    pub max_force: f64,
    /// Largest strain gradient component per atom, in eV/Angstrom.
    pub cell_force: f64,
    /// Structure at these variables.
    pub structure: Structure,
}

impl Evaluation {
    /// Trajectory record for this evaluation.
    pub(crate) fn frame(
        &self,
        step: usize,
        point: Option<usize>,
        cadence: TrajectoryOptions,
    ) -> TrajectoryFrame {
        TrajectoryFrame {
            step,
            point,
            time_fs: 0.0,
            energy: self.energy,
            enthalpy: self.enthalpy,
            max_force: self.max_force,
            volume: self.structure.lattice().volume().abs(),
            lattice: *self.structure.lattice(),
            positions: cadence.positions.then(|| {
                self.structure
                    .sites()
                    .iter()
                    .map(xtalforge_core::Site::frac)
                    .collect()
            }),
        }
    }
}

/// The energy landscape the minimizers walk.
pub(crate) struct Problem<'a> {
    calculator: &'a dyn Calculator,
    dof: Dof,
    /// External pressure in eV/Angstrom^3.
    pressure: f64,
}

impl<'a> Problem<'a> {
    /// Prepares the variables and the pressure term.
    pub(crate) fn new(
        calculator: &'a dyn Calculator,
        structure: &Structure,
        options: &OptimizeOptions,
    ) -> Result<Self, Error> {
        options.validate()?;
        let dof = Dof::new(
            structure,
            options.mode,
            options.cell,
            &options.fixed_atoms,
            options.cell_freedom,
        )?;
        Ok(Self {
            calculator,
            dof,
            pressure: options.pressure_gpa * GPA_TO_EV_PER_A3,
        })
    }

    /// Initial variable vector.
    pub(crate) fn start(&self) -> Vec<f64> {
        self.dof.start()
    }

    /// Projects a step direction onto the allowed subspace.
    pub(crate) fn project(&self, direction: &mut [f64]) {
        self.dof.project_direction(direction);
    }

    /// Energy, gradient and derived numbers at one variable vector.
    pub(crate) fn evaluate(&mut self, variables: &[f64]) -> Result<Evaluation, Error> {
        let structure = self.dof.structure_of(variables)?;
        let (energy, forces) = self.calculator.energy_forces(&structure)?;
        let volume = structure.lattice().volume().abs();
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
            enthalpy: energy + self.pressure * volume,
            gradient,
            max_force,
            cell_force,
            structure,
        })
    }

    /// Strains the cell from the calculator's stress, or by central differences
    /// when the engine provides none or the freedom is constrained.
    fn strain_gradient(
        &mut self,
        variables: &[f64],
        structure: &Structure,
    ) -> Result<[f64; 6], Error> {
        if self.dof.cell_freedom().is_full()
            && let Some(stress) = self.calculator.stress(structure)?
        {
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
            for (index, row) in converted.iter_mut().enumerate() {
                if let Some(value) = row.get_mut(index) {
                    *value += self.pressure;
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
            let energy_high = self.objective(&high)?;
            let energy_low = self.objective(&low)?;
            *slot = (energy_high - energy_low) / (2.0 * step);
        }
        Ok(gradient)
    }

    /// Objective at one variable vector, used by the difference gradient.
    fn objective(&mut self, variables: &[f64]) -> Result<f64, Error> {
        let structure = self.dof.structure_of(variables)?;
        let energy = self.calculator.energy(&structure)?;
        Ok(energy + self.pressure * structure.lattice().volume().abs())
    }

    /// Whether the force and energy criteria are met.
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
            .map(|energy| (evaluation.enthalpy - energy).abs() <= options.energy_tolerance)
            .unwrap_or(false)
    }

    /// Caps the largest per-site Cartesian displacement in a direction.
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

    /// Wraps one evaluation into a finished optimization result.
    pub(crate) fn result(
        &self,
        evaluation: Evaluation,
        steps: usize,
        converged: bool,
        cancelled: bool,
        trajectory: Trajectory,
    ) -> OptimizationResult {
        OptimizationResult {
            structure: evaluation.structure,
            energy: evaluation.energy,
            enthalpy: evaluation.enthalpy,
            max_force: evaluation.max_force,
            steps,
            converged,
            cancelled,
            trajectory,
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
