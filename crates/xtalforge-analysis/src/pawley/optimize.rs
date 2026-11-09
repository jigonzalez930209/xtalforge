//! Nelder-Mead search for the Pawley non-linear parameters (Task 7.6.2).

use xtalforge_core::Error;

pub(super) fn nelder_mead(
    start: &[f64],
    steps: &[f64],
    mut objective: impl FnMut(&[f64]) -> Result<f64, Error>,
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
    for _ in 0..200 {
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
        let reflected = centroid
            .iter()
            .zip(worst.iter())
            .map(|(center, value)| center + (center - value))
            .collect::<Vec<_>>();
        let reflected_value = objective(&reflected)?;
        let best_value = simplex.first().map_or(f64::INFINITY, |value| value.1);
        let second_worst = simplex
            .get(count.saturating_sub(1))
            .map_or(f64::INFINITY, |value| value.1);
        if reflected_value < best_value {
            let expanded = centroid
                .iter()
                .zip(worst.iter())
                .map(|(center, value)| center + 2.0 * (center - value))
                .collect::<Vec<_>>();
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
            let contracted = centroid
                .iter()
                .zip(worst.iter())
                .map(|(center, value)| center + 0.5 * (value - center))
                .collect::<Vec<_>>();
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
                    let shrunk = anchor
                        .iter()
                        .zip(point.iter())
                        .map(|(base, value)| base + 0.5 * (value - base))
                        .collect::<Vec<_>>();
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
