//! Active-learning selection strategies (Task 21.4.1).

#![allow(clippy::indexing_slicing)] // bounded candidate indices

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

/// Selection strategy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SelectionStrategy {
    /// Highest predictive uncertainty.
    Uncertainty,
    /// Farthest point in descriptor space.
    FarthestPoint,
    /// Uncertainty first, then diversity within the batch.
    Diversity,
}

/// Options for [`select_batch`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SelectionOptions {
    /// Number of candidates to select.
    pub batch: usize,
    /// Strategy.
    pub strategy: SelectionStrategy,
    /// Minimum Euclidean distance for the diversity strategy; zero disables.
    pub diversity_radius: f64,
}

impl Default for SelectionOptions {
    fn default() -> Self {
        Self {
            batch: 8,
            strategy: SelectionStrategy::Uncertainty,
            diversity_radius: 0.0,
        }
    }
}

/// Selects a batch of candidate indices.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a zero batch, mismatched lengths or
/// non-finite values.
pub fn select_batch(
    features: &[Vec<f64>],
    uncertainties: &[f64],
    options: &SelectionOptions,
) -> Result<Vec<usize>, Error> {
    if options.batch == 0 {
        return Err(Error::InvalidInput(
            "the selection batch must be at least one".to_owned(),
        ));
    }
    if features.len() != uncertainties.len() {
        return Err(Error::InvalidInput(
            "features and uncertainties must have the same length".to_owned(),
        ));
    }
    if !uncertainties.iter().all(|value| value.is_finite()) {
        return Err(Error::InvalidInput(
            "uncertainties must be finite".to_owned(),
        ));
    }
    if features
        .iter()
        .any(|row| row.iter().any(|value| !value.is_finite()))
    {
        return Err(Error::InvalidInput("features must be finite".to_owned()));
    }
    let n = features.len();
    if n == 0 {
        return Ok(Vec::new());
    }
    let batch = options.batch.min(n);
    let mut selected = match options.strategy {
        SelectionStrategy::Uncertainty => {
            let mut order: Vec<usize> = (0..n).collect();
            order.sort_by(|left, right| {
                uncertainties[*right]
                    .total_cmp(&uncertainties[*left])
                    .then_with(|| left.cmp(right))
            });
            order.truncate(batch);
            order
        }
        SelectionStrategy::FarthestPoint => farthest_point(features, batch),
        SelectionStrategy::Diversity => {
            let mut order: Vec<usize> = (0..n).collect();
            order.sort_by(|left, right| {
                uncertainties[*right]
                    .total_cmp(&uncertainties[*left])
                    .then_with(|| left.cmp(right))
            });
            let mut selected: Vec<usize> = Vec::new();
            for index in order {
                if selected.len() >= batch {
                    break;
                }
                let candidate: &[f64] = features[index].as_slice();
                let far_enough = selected.iter().all(|other| {
                    let representative: &[f64] = features[*other].as_slice();
                    distance(candidate, representative) >= options.diversity_radius
                });
                if far_enough || selected.is_empty() {
                    selected.push(index);
                }
            }
            selected
        }
    };
    selected.sort_unstable();
    Ok(selected)
}

/// Greedy farthest-point selection, seeded by the highest uncertainty.
#[must_use]
pub fn farthest_point(features: &[Vec<f64>], batch: usize) -> Vec<usize> {
    if features.is_empty() || batch == 0 {
        return Vec::new();
    }
    let mut selected = vec![0usize];
    while selected.len() < batch.min(features.len()) {
        let mut best = None;
        let mut best_distance = -1.0;
        for index in 0..features.len() {
            if selected.contains(&index) {
                continue;
            }
            let nearest = selected
                .iter()
                .map(|other| distance(features[index].as_slice(), features[*other].as_slice()))
                .fold(f64::INFINITY, f64::min);
            if nearest > best_distance {
                best_distance = nearest;
                best = Some(index);
            }
        }
        let Some(index) = best else {
            break;
        };
        selected.push(index);
    }
    selected.sort_unstable();
    selected
}

/// Sorts candidates by descending uncertainty.
#[must_use]
pub fn batch_uncertainty(uncertainties: &[f64]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..uncertainties.len()).collect();
    order.sort_by(|left, right| {
        uncertainties[*right]
            .total_cmp(&uncertainties[*left])
            .then_with(|| left.cmp(right))
    });
    order
}

fn distance(first: &[f64], second: &[f64]) -> f64 {
    first
        .iter()
        .zip(second)
        .map(|(a, b)| (a - b) * (a - b))
        .sum::<f64>()
        .sqrt()
}
