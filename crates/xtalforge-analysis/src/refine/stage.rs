//! One bounded refinement stage of the Rietveld fit (Task 7.8.1).
//!
//! Every stage uses a Nelder-Mead simplex on the stage's free parameters with
//! the scale and background refitted analytically at each vertex. A direct
//! search is used instead of Gauss-Newton because the zero shift, the cell and
//! the width are strongly correlated with each other and with the scale: the
//! normal-equation direction at a displaced start can point uphill, while the
//! simplex still finds the valley. The stages keep the McCusker order; the
//! document in [`super`] records the deviations.

#![allow(clippy::indexing_slicing)] // bounded local simplex

use super::RefineOptions;
use super::model::{Model, ModelContext, bounds, clamp};
use xtalforge_core::Error;

/// Result of one stage.
pub(super) struct StageResult {
    /// Final model of the stage.
    pub(super) model: Model,
    /// Whether the stage reached a stationary point.
    pub(super) converged: bool,
    /// Simplex iterations used by the stage.
    pub(super) iterations: usize,
}

/// Refines the parameters at `indices` while every other value stays fixed.
pub(super) fn run_stage(
    context: &ModelContext,
    values: &mut [f64],
    terms: usize,
    indices: &[usize],
    options: &RefineOptions,
) -> Result<StageResult, Error> {
    let mut current = context.model(values, terms)?;
    if indices.is_empty() {
        return Ok(StageResult {
            model: current,
            converged: true,
            iterations: 0,
        });
    }
    let start = indices
        .iter()
        .filter_map(|index| values.get(*index).copied())
        .collect::<Vec<_>>();
    let steps = indices
        .iter()
        .map(|index| {
            context
                .parameters
                .get(*index)
                .map_or(0.01, |parameter| step(parameter.kind))
        })
        .collect::<Vec<_>>();
    let objective = |free: &[f64]| -> Result<f64, Error> {
        let mut trial = values.to_vec();
        for (column, index) in indices.iter().enumerate() {
            let Some(parameter) = context.parameters.get(*index) else {
                continue;
            };
            let (low, high) = bounds(parameter.kind);
            let value = clamp(free.get(column).copied().unwrap_or(0.0), low, high);
            if let Some(slot) = trial.get_mut(*index) {
                *slot = value;
            }
        }
        Ok(context.model(&trial, terms)?.rwp)
    };
    let iterations = options.max_iterations.saturating_mul(4).clamp(40, 400);
    let (best, _, converged) = nelder_mead(&start, &steps, objective, iterations)?;
    for (column, index) in indices.iter().enumerate() {
        let Some(parameter) = context.parameters.get(*index) else {
            continue;
        };
        let (low, high) = bounds(parameter.kind);
        if let Some(slot) = values.get_mut(*index) {
            *slot = clamp(best.get(column).copied().unwrap_or(*slot), low, high);
        }
    }
    current = context.model(values, terms)?;
    Ok(StageResult {
        model: current,
        converged,
        iterations,
    })
}

/// Starting simplex edge for one parameter kind.
fn step(kind: super::model::ParameterKind) -> f64 {
    use super::model::ParameterKind;
    match kind {
        ParameterKind::Position { .. } => 0.005,
        ParameterKind::AdpScale => 0.1,
        ParameterKind::Occupancy { .. } => 0.03,
        ParameterKind::Zero => 0.01,
        ParameterKind::Cell { index } if index < 3 => 0.02,
        ParameterKind::Cell { .. } => 0.05,
        ParameterKind::WidthScale | ParameterKind::Texture => 0.03,
    }
}

/// Nelder-Mead search over the free values of one stage.
fn nelder_mead(
    start: &[f64],
    steps: &[f64],
    mut objective: impl FnMut(&[f64]) -> Result<f64, Error>,
    iterations: usize,
) -> Result<(Vec<f64>, f64, bool), Error> {
    let count = start.len();
    let mut simplex: Vec<(Vec<f64>, f64)> = Vec::with_capacity(count + 1);
    simplex.push((start.to_vec(), objective(start)?));
    for column in 0..count {
        let mut point = start.to_vec();
        if let Some(slot) = point.get_mut(column) {
            *slot += steps.get(column).copied().unwrap_or(0.01);
        }
        let value = objective(&point)?;
        simplex.push((point, value));
    }
    let mut converged = false;
    for _ in 0..iterations {
        simplex.sort_by(|first, second| first.1.total_cmp(&second.1));
        let best = simplex.first().map_or(0.0, |value| value.1);
        let spread = simplex
            .iter()
            .map(|value| (value.1 - best).abs())
            .fold(0.0_f64, f64::max);
        if spread <= 1e-10 * best.abs().max(1.0) {
            converged = true;
            break;
        }
        let centroid: Vec<f64> = (0..count)
            .map(|column| {
                simplex
                    .iter()
                    .take(count)
                    .map(|value| value.0.get(column).copied().unwrap_or(0.0))
                    .sum::<f64>()
                    / count as f64
            })
            .collect();
        let worst = simplex
            .get(count)
            .map(|value| value.0.clone())
            .unwrap_or_default();
        let reflected: Vec<f64> = centroid
            .iter()
            .zip(worst.iter())
            .map(|(center, value)| center + (center - value))
            .collect();
        let reflected_value = objective(&reflected)?;
        let best_value = simplex.first().map_or(f64::INFINITY, |value| value.1);
        let second_worst = simplex
            .get(count.saturating_sub(1))
            .map_or(f64::INFINITY, |value| value.1);
        if reflected_value < best_value {
            let expanded: Vec<f64> = centroid
                .iter()
                .zip(worst.iter())
                .map(|(center, value)| center + 2.0 * (center - value))
                .collect();
            let expanded_value = objective(&expanded)?;
            if expanded_value < reflected_value {
                simplex[count] = (expanded, expanded_value);
            } else {
                simplex[count] = (reflected, reflected_value);
            }
        } else if reflected_value < second_worst {
            simplex[count] = (reflected, reflected_value);
        } else {
            let worst_value = simplex.get(count).map_or(f64::INFINITY, |value| value.1);
            let contracted: Vec<f64> = centroid
                .iter()
                .zip(worst.iter())
                .map(|(center, value)| center + 0.5 * (value - center))
                .collect();
            let contracted_value = objective(&contracted)?;
            if contracted_value < worst_value {
                simplex[count] = (contracted, contracted_value);
            } else {
                let anchor = simplex
                    .first()
                    .map(|value| value.0.clone())
                    .unwrap_or_default();
                for index in 1..=count {
                    let Some((point, _)) = simplex.get(index) else {
                        continue;
                    };
                    let shrunk: Vec<f64> = anchor
                        .iter()
                        .zip(point.iter())
                        .map(|(base, value)| base + 0.5 * (value - base))
                        .collect();
                    let value = objective(&shrunk)?;
                    simplex[index] = (shrunk, value);
                }
            }
        }
    }
    simplex.sort_by(|first, second| first.1.total_cmp(&second.1));
    let (point, value) = simplex.into_iter().next().unwrap_or_default();
    Ok((point, value, converged))
}
