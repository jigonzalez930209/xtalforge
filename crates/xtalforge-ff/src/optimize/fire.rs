//! FIRE minimizer (Task 4.4.3).
//!
//! Bitzek et al., *Phys. Rev. Lett.* **2006**, 97, 170201.  Molecular dynamics
//! with the velocity steered towards the force and the time step grown while
//! progress continues, thrown away the moment it stops.  Slower than L-BFGS
//! near a minimum and far more forgiving a long way from one, which is what a
//! hand-built structure needs.
//!
//! [`fire_with`] adds cancellation, progress and a trajectory; [`fire`] is the
//! same run without a control object.

use xtalforge_core::Structure;

use crate::Error;
use crate::calculator::Calculator;
use crate::optimize::{OptimizationResult, OptimizeOptions, Problem, axpy, dot, norm, record};
use crate::run::{ProgressEvent, RunControl};
use crate::trajectory::Trajectory;

/// Minimizes a structure with FIRE.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for bad options or variable setup, and any
/// typing, topology or data error from the calculator.
pub fn fire(
    calculator: &dyn Calculator,
    structure: &Structure,
    options: &OptimizeOptions,
) -> Result<OptimizationResult, Error> {
    fire_with(calculator, structure, options, &mut RunControl::new())
}

/// Minimizes a structure with FIRE, reporting progress and honouring
/// cancellation.
///
/// A cancelled run returns the last geometry with `converged = false` and
/// `cancelled = true`; it is not an error.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for bad options or variable setup, and any
/// typing, topology or data error from the calculator.
pub fn fire_with(
    calculator: &dyn Calculator,
    structure: &Structure,
    options: &OptimizeOptions,
    control: &mut RunControl<'_>,
) -> Result<OptimizationResult, Error> {
    options.validate()?;
    let mut problem = Problem::new(calculator, structure, options)?;
    let mut variables = problem.start();
    let mut evaluation = problem.evaluate(&variables)?;
    let mut trajectory = Trajectory::new();
    record(&mut trajectory, options.trajectory, 0, None, &evaluation);
    control.report(ProgressEvent::Optimizer {
        step: 0,
        energy: evaluation.energy,
        max_force: evaluation.max_force,
    });
    if control.is_cancelled() {
        return Ok(problem.result(evaluation, 0, false, true, trajectory));
    }
    if problem.converged(&evaluation, None, options) {
        return Ok(problem.result(evaluation, 0, true, false, trajectory));
    }

    let mut velocity = vec![0.0; variables.len()];
    let mut alpha = 0.1;
    let mut dt: f64 = 0.1;
    let dt_max = 1.0;
    let n_min = 5;
    let f_inc = 1.1;
    let f_dec = 0.5;
    let alpha_start = 0.1;
    let f_alpha = 0.99;
    let mut since_uphill = 0_usize;

    for step in 1..=options.max_steps {
        if control.is_cancelled() {
            return Ok(problem.result(evaluation, step - 1, false, true, trajectory));
        }
        let force: Vec<f64> = evaluation.gradient.iter().map(|value| -value).collect();
        let power = dot(&force, &velocity);
        if power > 0.0 {
            let force_norm = norm(&force);
            let velocity_norm = norm(&velocity);
            if force_norm > 0.0 {
                for (index, value) in velocity.iter_mut().enumerate() {
                    let target = force.get(index).copied().unwrap_or(0.0);
                    *value = (1.0 - alpha) * *value + alpha * velocity_norm * target / force_norm;
                }
            }
            if since_uphill > n_min {
                dt = (dt * f_inc).min(dt_max);
                alpha *= f_alpha;
            }
            since_uphill += 1;
        } else {
            velocity.fill(0.0);
            since_uphill = 0;
            dt *= f_dec;
            alpha = alpha_start;
        }

        for (value, component) in velocity.iter_mut().zip(&force) {
            *value += dt * component;
        }
        let mut direction: Vec<f64> = velocity.iter().map(|value| dt * value).collect();
        problem.cap_step(&mut direction, options.max_step);
        problem.project(&mut direction);

        let mut trial = variables.clone();
        axpy(&mut trial, 1.0, &direction);
        let new_evaluation = problem.evaluate(&trial)?;
        record(
            &mut trajectory,
            options.trajectory,
            step,
            None,
            &new_evaluation,
        );
        control.report(ProgressEvent::Optimizer {
            step,
            energy: new_evaluation.energy,
            max_force: new_evaluation.max_force,
        });
        let converged = problem.converged(&new_evaluation, Some(evaluation.enthalpy), options);
        if converged {
            return Ok(problem.result(new_evaluation, step, true, false, trajectory));
        }

        variables = trial;
        evaluation = new_evaluation;
        // The direction is `dt * velocity`; carrying it back as a velocity
        // keeps the next acceleration in the same units as this one.
        velocity = direction.iter().map(|value| value / dt).collect();
    }

    Ok(problem.result(evaluation, options.max_steps, false, false, trajectory))
}
