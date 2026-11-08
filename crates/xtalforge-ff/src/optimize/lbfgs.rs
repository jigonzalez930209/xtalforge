//! L-BFGS minimizer with an Armijo backtracking line search (Task 4.4.3).
//!
//! The limited-memory quasi-Newton method of Liu and Nocedal: the last ten
//! curvature pairs stand in for the inverse Hessian and the two-loop recursion
//! applies them.  Pairs that fail `s . y > 0` are dropped rather than stored,
//! because keeping one makes the approximation indefinite and the next
//! direction points uphill.  A failed line search drops the history and tries
//! steepest descent once from the same point; only then does the run stop, and
//! it says it did not converge.

use xtalforge_core::Structure;

use crate::Error;
use crate::calculator::Calculator;
use crate::optimize::{
    Evaluation, OptimizationResult, OptimizeOptions, Problem, axpy, difference, dot,
};

/// Curvature pairs kept between restarts.
const MEMORY: usize = 10;

/// Minimizes a structure with L-BFGS.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for bad options or variable setup, and any
/// typing, topology or data error from the calculator.
pub fn lbfgs(
    calculator: &dyn Calculator,
    structure: &Structure,
    options: &OptimizeOptions,
) -> Result<OptimizationResult, Error> {
    options.validate()?;
    let mut problem = Problem::new(calculator, structure, options)?;
    let mut variables = problem.start();
    let mut evaluation = problem.evaluate(&variables)?;
    if problem.converged(&evaluation, None, options) {
        return Ok(problem.result(evaluation, 0, true));
    }

    let mut history: Vec<(Vec<f64>, Vec<f64>, f64)> = Vec::new();
    for step in 1..=options.max_steps {
        let mut direction = two_loop(&evaluation.gradient, &history);
        problem.project(&mut direction);
        if dot(&direction, &evaluation.gradient) >= 0.0 {
            history.clear();
            direction = evaluation.gradient.iter().map(|value| -value).collect();
            problem.project(&mut direction);
        }

        let mut accepted = line_search(
            &mut problem,
            &variables,
            &evaluation,
            &direction,
            options.max_step,
        )?;
        if accepted.is_none() && !history.is_empty() {
            history.clear();
            let mut steepest: Vec<f64> = evaluation.gradient.iter().map(|value| -value).collect();
            problem.project(&mut steepest);
            accepted = line_search(
                &mut problem,
                &variables,
                &evaluation,
                &steepest,
                options.max_step,
            )?;
        }
        let Some((new_variables, new_evaluation)) = accepted else {
            return Ok(problem.result(evaluation, step, false));
        };

        let s = difference(&new_variables, &variables);
        let y = difference(&new_evaluation.gradient, &evaluation.gradient);
        let curvature = dot(&s, &y);
        if curvature > 1e-12 {
            history.push((s, y, 1.0 / curvature));
            if history.len() > MEMORY {
                history.remove(0);
            }
        }

        let converged = problem.converged(&new_evaluation, Some(evaluation.energy), options);
        variables = new_variables;
        evaluation = new_evaluation;
        if converged {
            return Ok(problem.result(evaluation, step, true));
        }
    }

    Ok(problem.result(evaluation, options.max_steps, false))
}

/// The two-loop recursion: `H * gradient` from the stored curvature pairs.
fn two_loop(gradient: &[f64], history: &[(Vec<f64>, Vec<f64>, f64)]) -> Vec<f64> {
    let mut q = gradient.to_vec();
    let mut coefficients = Vec::with_capacity(history.len());
    for (s, y, rho) in history.iter().rev() {
        let coefficient = rho * dot(s, &q);
        axpy(&mut q, -coefficient, y);
        coefficients.push(coefficient);
    }
    let mut scale = 1.0;
    if let Some((s, y, _)) = history.last() {
        let denominator = dot(y, y);
        if denominator > 0.0 {
            scale = dot(s, y) / denominator;
        }
    }
    if !scale.is_finite() || scale <= 0.0 {
        scale = 1.0;
    }
    for value in &mut q {
        *value *= scale;
    }
    for ((s, y, rho), coefficient) in history.iter().zip(coefficients.iter().rev()) {
        let beta = rho * dot(y, &q);
        axpy(&mut q, coefficient - beta, s);
    }
    q
}

/// Backtracks along `direction` until the energy has gone down enough.
fn line_search(
    problem: &mut Problem<'_>,
    variables: &[f64],
    evaluation: &Evaluation,
    direction: &[f64],
    max_step: f64,
) -> Result<Option<(Vec<f64>, Evaluation)>, Error> {
    let slope = dot(&evaluation.gradient, direction);
    if slope >= 0.0 {
        return Ok(None);
    }
    let mut length = 1.0;
    let longest = longest_site_step(direction);
    if longest * length > max_step && longest > 0.0 {
        length = max_step / longest;
    }
    for _attempt in 0..20 {
        let mut trial = variables.to_vec();
        axpy(&mut trial, length, direction);
        let new_evaluation = problem.evaluate(&trial)?;
        if new_evaluation.energy <= evaluation.energy + 1e-4 * length * slope {
            return Ok(Some((trial, new_evaluation)));
        }
        let curvature =
            (new_evaluation.energy - evaluation.energy - slope * length) / (length * length);
        let guess = if curvature > 0.0 {
            -slope / (2.0 * curvature)
        } else {
            0.5 * length
        };
        length = guess.clamp(0.1 * length, 0.5 * length);
    }
    Ok(None)
}

/// Largest per-site Cartesian displacement in a flat direction.
fn longest_site_step(direction: &[f64]) -> f64 {
    let mut longest: f64 = 0.0;
    let mut index = 0;
    while index + 2 < direction.len() {
        let vector = [
            direction.get(index).copied().unwrap_or(0.0),
            direction.get(index + 1).copied().unwrap_or(0.0),
            direction.get(index + 2).copied().unwrap_or(0.0),
        ];
        longest = longest.max(crate::math::norm(vector));
        index += 3;
    }
    longest
}
