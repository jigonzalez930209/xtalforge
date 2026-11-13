//! Pareto front, sweeps and knee detection for combined objectives
//! (Tasks 7.9.2 and 14.1.5).
//!
//! [`pareto_front`] is the original non-dominated filter. [`pareto_sweep`]
//! adds the sweep that a weight scan needs: evaluate a trade-off curve, mark
//! the non-dominated points and name the knee, the point of the front farthest
//! from the chord joining its two ends.

use xtalforge_core::Error;

/// Indices of the non-dominated points, minimising both coordinates.
///
/// A point is dominated when another point is no worse in both objectives and
/// strictly better in at least one. The result is sorted by the first
/// objective and then the second.
#[must_use]
pub fn pareto_front(points: &[(f64, f64)]) -> Vec<usize> {
    let mut front = Vec::new();
    for (index, (first, second)) in points.iter().enumerate() {
        let dominated = points.iter().enumerate().any(|(other, (x, y))| {
            other != index && x <= first && y <= second && (x < first || y < second)
        });
        if !dominated {
            front.push(index);
        }
    }
    front.sort_by(|a, b| {
        let first = points
            .get(*a)
            .copied()
            .unwrap_or((f64::INFINITY, f64::INFINITY));
        let second = points
            .get(*b)
            .copied()
            .unwrap_or((f64::INFINITY, f64::INFINITY));
        first
            .0
            .partial_cmp(&second.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                first
                    .1
                    .partial_cmp(&second.1)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });
    front
}

/// One sample of a weight sweep.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ParetoSample {
    /// Weight of the first objective in `[0, 1]`.
    pub weight: f64,
    /// First objective (for example the energy term).
    pub first: f64,
    /// Second objective (for example Rwp).
    pub second: f64,
    /// Linear scalarization `weight * first + (1 - weight) * second`.
    pub combined: f64,
}

/// A weight sweep with its front and knee.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct ParetoSweep {
    /// Every evaluated sample, in weight order.
    pub samples: Vec<ParetoSample>,
    /// Indices of the non-dominated samples, sorted by the first objective.
    pub front: Vec<usize>,
    /// Index of the knee among `samples`, when the front has an interior bend.
    pub knee: Option<usize>,
}

/// Evaluates a trade-off curve over weights, marks the front and finds the
/// knee.
///
/// `evaluate` returns `(first, second)` for one weight; both values must be
/// finite. The combined value is the linear scalarization with the first
/// objective weighted by `weight`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty or non-finite weight list or a
/// non-finite objective returned by `evaluate`.
pub fn pareto_sweep<F>(weights: &[f64], mut evaluate: F) -> Result<ParetoSweep, Error>
where
    F: FnMut(f64) -> Result<(f64, f64), Error>,
{
    if weights.is_empty() {
        return Err(Error::InvalidInput(
            "a Pareto sweep needs at least one weight".to_owned(),
        ));
    }
    let mut samples = Vec::with_capacity(weights.len());
    for weight in weights {
        if !weight.is_finite() || !(0.0..=1.0).contains(weight) {
            return Err(Error::InvalidInput(format!(
                "sweep weight must be in [0, 1], got {weight}"
            )));
        }
        let (first, second) = evaluate(*weight)?;
        if !first.is_finite() || !second.is_finite() {
            return Err(Error::InvalidInput(format!(
                "objective at weight {weight} is not finite: ({first}, {second})"
            )));
        }
        samples.push(ParetoSample {
            weight: *weight,
            first,
            second,
            combined: weight * first + (1.0 - weight) * second,
        });
    }
    let points: Vec<(f64, f64)> = samples
        .iter()
        .map(|sample| (sample.first, sample.second))
        .collect();
    let front = pareto_front(&points);
    let knee = knee_index(
        &front
            .iter()
            .filter_map(|index| points.get(*index).copied())
            .collect::<Vec<_>>(),
    )
    .and_then(|position| front.get(position).copied());
    Ok(ParetoSweep {
        samples,
        front,
        knee,
    })
}

/// Index of the knee of a Pareto front, or `None` when there is no bend.
///
/// The points are expected in front order (increasing first objective). Both
/// axes are normalised to `[0, 1]`; the knee is the point with the largest
/// perpendicular distance from the chord joining the first and last point. A
/// straight front has no knee and returns `None`.
#[must_use]
pub fn knee_index(points: &[(f64, f64)]) -> Option<usize> {
    if points.len() < 3 {
        return None;
    }
    let first_min = points
        .iter()
        .map(|point| point.0)
        .fold(f64::INFINITY, f64::min);
    let first_max = points
        .iter()
        .map(|point| point.0)
        .fold(f64::NEG_INFINITY, f64::max);
    let second_min = points
        .iter()
        .map(|point| point.1)
        .fold(f64::INFINITY, f64::min);
    let second_max = points
        .iter()
        .map(|point| point.1)
        .fold(f64::NEG_INFINITY, f64::max);
    let first_span = first_max - first_min;
    let second_span = second_max - second_min;
    if first_span <= 0.0 || second_span <= 0.0 {
        return None;
    }
    let normalised: Vec<(f64, f64)> = points
        .iter()
        .map(|point| {
            (
                (point.0 - first_min) / first_span,
                (point.1 - second_min) / second_span,
            )
        })
        .collect();
    let (first, second) = (normalised.first().copied()?, normalised.last().copied()?);
    let (dx, dy) = (second.0 - first.0, second.1 - first.1);
    let length = (dx * dx + dy * dy).sqrt();
    if length <= 1e-12 {
        return None;
    }
    let mut best = None;
    let mut best_distance = 0.0;
    for (index, point) in normalised.iter().enumerate() {
        let distance = ((dy * (point.0 - first.0) - dx * (point.1 - first.1)) / length).abs();
        if distance > best_distance {
            best_distance = distance;
            best = Some(index);
        }
    }
    if best_distance <= 1e-12 { None } else { best }
}
