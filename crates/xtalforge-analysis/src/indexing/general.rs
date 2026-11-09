//! Metric-family search behind [`super::index_pattern`] (Task 7.5.1).

#![allow(clippy::indexing_slicing)] // fixed-size local matrices

use std::collections::BTreeSet;

use xtalforge_core::Error;

use crate::peaks::Peak;

use super::fom::{evaluate_cell, index_limit, peak_q, two_theta_of_q};
use super::{CellMetric, IndexedCell, IndexingOptions, LatticeSystem};

/// `Q = A h^2 + B k^2 + C l^2 + D h l` coefficients and the Miller indices
/// they came from.
type Combo = ([f64; 4], [i32; 3]);

/// Solves every requested family and returns the raw candidates.
pub(crate) fn index_general(
    peaks: &[Peak],
    wavelength: f64,
    options: &IndexingOptions,
) -> Result<Vec<IndexedCell>, Error> {
    let mut found = Vec::new();
    for system in &options.systems {
        let mut seen = BTreeSet::new();
        let solved = solve_family(*system, peaks, wavelength, options, &mut seen)?;
        if solved.is_empty() && system.parameters() >= 3 {
            // A larger first-index window is only paid for when the cheap one
            // found nothing.
            let wider = IndexingOptions {
                max_index: options.max_index + 1,
                ..options.clone()
            };
            let mut retry_seen = BTreeSet::new();
            found.extend(solve_family(
                *system,
                peaks,
                wavelength,
                &wider,
                &mut retry_seen,
            )?);
        } else {
            found.extend(solved);
        }
    }
    Ok(found)
}

fn solve_family(
    system: LatticeSystem,
    peaks: &[Peak],
    wavelength: f64,
    options: &IndexingOptions,
    seen: &mut BTreeSet<[i64; 6]>,
) -> Result<Vec<IndexedCell>, Error> {
    let parameters = system.parameters();
    if peaks.len() < parameters {
        return Ok(Vec::new());
    }
    let combos = combinations(system, effective_index(system, options.max_index));
    if combos.is_empty() {
        return Ok(Vec::new());
    }
    // A monoclinic solve needs a reflection with both h and l non-zero to
    // reach the cross term, so four of the first six lines are tried; the
    // first four alone are linearly dependent.
    // The first lines of a family are often degenerate for the solve (a
    // tetragonal pattern starts with h00 and hk0, both with l = 0), so a few
    // extra lines are offered and the singular subsets are skipped.
    let window = match system {
        LatticeSystem::Cubic => 3.min(peaks.len()),
        LatticeSystem::Monoclinic => 6.min(peaks.len()),
        _ => 4.min(peaks.len()),
    };
    let mut found = Vec::new();
    let mut assignment = vec![0usize; parameters];
    let subsets = choose(window, parameters);
    let per_subset = (options.max_evaluations / subsets.len().max(1)).max(1);
    for subset in subsets {
        let mut budget = per_subset;
        assign_combos(
            0,
            &subset,
            &combos,
            peaks,
            wavelength,
            options,
            system,
            seen,
            &mut assignment,
            &mut budget,
            &mut found,
        );
    }
    Ok(found)
}

fn effective_index(system: LatticeSystem, max_index: i32) -> i32 {
    match system {
        // Four unknowns make the enumeration the expensive part; the wide
        // retry (options.max_index + 1) is the one that reaches index two.
        LatticeSystem::Monoclinic if max_index <= 3 => 1,
        LatticeSystem::Monoclinic => 2,
        _ => max_index,
    }
}

#[allow(clippy::too_many_arguments)] // the recursion threads one search context
fn assign_combos(
    depth: usize,
    subset: &[usize],
    combos: &[Combo],
    peaks: &[Peak],
    wavelength: f64,
    options: &IndexingOptions,
    system: LatticeSystem,
    seen: &mut BTreeSet<[i64; 6]>,
    assignment: &mut [usize],
    budget: &mut usize,
    found: &mut Vec<IndexedCell>,
) {
    if *budget == 0 {
        return;
    }
    if depth == assignment.len() {
        *budget -= 1;
        evaluate(
            system, subset, combos, assignment, peaks, wavelength, options, seen, found,
        );
        return;
    }
    for index in 0..combos.len() {
        if let Some(slot) = assignment.get_mut(depth) {
            *slot = index;
        }
        assign_combos(
            depth + 1,
            subset,
            combos,
            peaks,
            wavelength,
            options,
            system,
            seen,
            assignment,
            budget,
            found,
        );
    }
}

#[allow(clippy::too_many_arguments)] // one candidate carries its search context
fn evaluate(
    system: LatticeSystem,
    subset: &[usize],
    combos: &[Combo],
    assignment: &[usize],
    peaks: &[Peak],
    wavelength: f64,
    options: &IndexingOptions,
    seen: &mut BTreeSet<[i64; 6]>,
    found: &mut Vec<IndexedCell>,
) {
    let parameters = system.parameters();
    let mut matrix = [[0.0; 4]; 4];
    let mut rhs = [0.0; 4];
    for (row, (combo_index, peak_index)) in assignment.iter().zip(subset).enumerate() {
        let Some(combo) = combos.get(*combo_index) else {
            return;
        };
        let Some(peak) = peaks.get(*peak_index) else {
            return;
        };
        let Ok(q) = peak_q(peak, wavelength) else {
            return;
        };
        if let Some(target) = matrix.get_mut(row) {
            *target = combo.0;
        }
        if let Some(slot) = rhs.get_mut(row) {
            *slot = q;
        }
    }
    let Some(solution) = solve(&matrix, &rhs, parameters) else {
        return;
    };
    let Some(cell) = system.cell(&solution) else {
        return;
    };
    if !cell.is_plausible(options.max_volume) {
        return;
    }
    let key = cell.key();
    if !seen.insert(key) {
        return;
    }
    let Some(next) = peaks
        .iter()
        .enumerate()
        .find(|(index, _)| !subset.contains(index))
    else {
        return;
    };
    if !predicts(&cell, next.1, wavelength, options.tolerance) {
        return;
    }
    let Ok(mut evaluation) = evaluate_cell(&cell, peaks, wavelength, options.tolerance) else {
        return;
    };
    let mut cell = cell;
    // Linear least-squares passes on the matched reflections remove the error
    // carried by the solving lines. The window grows from the unambiguous
    // low-angle lines to the whole pattern, and every pass re-matches the hkl;
    // that is what keeps a large cell from locking onto a wrong assignment at
    // high angle.
    for window in [8_usize, 16, usize::MAX] {
        let take = evaluation.indexed.len().min(window);
        let Some(refined) = refine(system, &evaluation.indexed[..take], wavelength) else {
            break;
        };
        if !refined.is_plausible(options.max_volume) {
            break;
        }
        let Ok(refined_evaluation) = evaluate_cell(&refined, peaks, wavelength, options.tolerance)
        else {
            break;
        };
        if refined_evaluation.completeness < evaluation.completeness {
            break;
        }
        cell = refined;
        evaluation = refined_evaluation;
    }
    if evaluation.completeness < options.min_completeness {
        return;
    }
    found.push(IndexedCell {
        system,
        cell,
        peaks: evaluation.indexed,
        completeness: evaluation.completeness,
        mean_error: evaluation.mean_error,
        m20: evaluation.figures.m20,
        f_n: evaluation.figures.f_n,
        coverage: evaluation.figures.coverage,
    });
}

/// Whether a cell predicts `peak` within `tolerance`.
///
/// The index range comes from the cell length and the peak's `Q`, so a large
/// cell is checked over a handful of indices instead of a fixed cube.
fn predicts(cell: &CellMetric, peak: &Peak, wavelength: f64, tolerance: f64) -> bool {
    let Ok(q) = peak_q(peak, wavelength) else {
        return false;
    };
    let max_index = index_limit(cell, q);
    for h in 0..=max_index {
        for k in 0..=max_index {
            for l in 0..=max_index {
                let hkl = [h, k, l];
                let Some(q) = cell.reciprocal_squared(hkl) else {
                    continue;
                };
                if (two_theta_of_q(q, wavelength) - peak.two_theta).abs() <= tolerance {
                    return true;
                }
            }
        }
    }
    false
}

/// Coefficients of `Q = A h^2 + B k^2 + C l^2 + D h l` for one family.
fn coefficients(system: LatticeSystem, hkl: [i32; 3]) -> [f64; 4] {
    let [h, k, l] = hkl.map(f64::from);
    match system {
        LatticeSystem::Cubic => [h * h + k * k + l * l, 0.0, 0.0, 0.0],
        LatticeSystem::Tetragonal => [h * h + k * k, l * l, 0.0, 0.0],
        LatticeSystem::Hexagonal => [h * h + h * k + k * k, l * l, 0.0, 0.0],
        LatticeSystem::Orthorhombic => [h * h, k * k, l * l, 0.0],
        LatticeSystem::Monoclinic => [h * h, k * k, l * l, h * l],
    }
}

/// Refines the metric by linear least squares on the matched reflections.
fn refine(
    system: LatticeSystem,
    indexed: &[super::IndexedPeak],
    wavelength: f64,
) -> Option<CellMetric> {
    let count = system.parameters();
    let mut normal = [[0.0; 4]; 4];
    let mut rhs = [0.0; 4];
    for peak in indexed {
        let theta = peak.two_theta.to_radians() / 2.0;
        let q = (2.0 * theta.sin() / wavelength).powi(2);
        let coefficients = coefficients(system, peak.hkl);
        for row in 0..count {
            for column in 0..count {
                normal[row][column] += coefficients[row] * coefficients[column];
            }
            rhs[row] += coefficients[row] * q;
        }
    }
    let solution = solve(&normal, &rhs, count)?;
    system.cell(&solution)
}

/// Distinct coefficient vectors for a family.
fn combinations(system: LatticeSystem, max_index: i32) -> Vec<Combo> {
    let mut seen: BTreeSet<[u64; 4]> = BTreeSet::new();
    let mut out: Vec<Combo> = Vec::new();
    let range: Vec<i32> = (-max_index..=max_index).collect();
    for h in &range {
        for k in 0..=max_index {
            for l in &range {
                if *h == 0 && k == 0 && *l == 0 {
                    continue;
                }
                let (h, k, l) = (*h, k, *l);
                let coefficients = coefficients(system, [h, k, l]);
                if seen.insert(coefficients.map(f64::to_bits)) {
                    out.push((coefficients, [h, k, l]));
                }
            }
        }
    }
    // Low-index assignments first, so a budget cut still tries the likely
    // cells before the exotic ones.
    out.sort_by_key(|(_, hkl)| hkl[0] * hkl[0] + hkl[1] * hkl[1] + hkl[2] * hkl[2]);
    out
}

/// Index subsets of `0..count` with `size` entries.
fn choose(count: usize, size: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    let mut current = Vec::with_capacity(size);
    choose_into(0, count, size, &mut current, &mut out);
    out
}

fn choose_into(
    start: usize,
    count: usize,
    size: usize,
    current: &mut Vec<usize>,
    out: &mut Vec<Vec<usize>>,
) {
    if current.len() == size {
        out.push(current.clone());
        return;
    }
    for index in start..count {
        current.push(index);
        choose_into(index + 1, count, size, current, out);
        current.pop();
    }
}

/// Solves the first `size` rows of an `n`-unknown linear system.
fn solve(matrix: &[[f64; 4]; 4], rhs: &[f64; 4], size: usize) -> Option<[f64; 4]> {
    let mut augmented = [[0.0; 5]; 4];
    for (row, target) in augmented.iter_mut().take(size).enumerate() {
        let source = matrix.get(row)?;
        for (column, cell) in target.iter_mut().take(size).enumerate() {
            *cell = *source.get(column)?;
        }
        target[size] = *rhs.get(row)?;
    }
    for pivot in 0..size {
        let mut best = pivot;
        for row in (pivot + 1)..size {
            if augmented[row][pivot].abs() > augmented[best][pivot].abs() {
                best = row;
            }
        }
        if augmented[best][pivot].abs() < 1e-12 {
            return None;
        }
        augmented.swap(pivot, best);
        let divisor = augmented[pivot][pivot];
        for cell in augmented[pivot].iter_mut().skip(pivot) {
            *cell /= divisor;
        }
        for row in 0..size {
            if row == pivot {
                continue;
            }
            let factor = augmented[row][pivot];
            if factor == 0.0 {
                continue;
            }
            let pivot_row = augmented[pivot];
            for (column, cell) in augmented[row].iter_mut().enumerate().skip(pivot) {
                *cell -= factor * pivot_row[column];
            }
        }
    }
    let mut solution = [0.0; 4];
    for row in 0..size {
        solution[row] = augmented[row][size];
    }
    Some(solution)
}
