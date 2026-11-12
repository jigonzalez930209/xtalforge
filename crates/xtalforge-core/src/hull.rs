#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Formation energies and convex-hull distances (Task 9.6.3).

use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};

use crate::Error;

/// One phase in a convex-hull calculation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PhaseEntry {
    /// Phase name.
    pub name: String,
    /// Element amounts in the formula unit.
    pub composition: BTreeMap<String, f64>,
    /// Total energy in eV for the formula unit.
    pub energy: f64,
}

/// One analysed phase.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PhasePoint {
    /// Phase name.
    pub name: String,
    /// Element fractions.
    pub composition: BTreeMap<String, f64>,
    /// Formation energy per atom in eV.
    pub formation_energy: f64,
    /// Energy of the lower convex hull at this composition, per atom.
    pub hull_energy: f64,
    /// Distance above the hull in eV/atom; zero for a stable phase.
    pub distance: f64,
    /// Whether the phase lies on the hull within `1e-6` eV/atom.
    pub stable: bool,
    /// Hull decomposition as `(entry index, weight)` for unstable phases.
    pub decomposition: Vec<(usize, f64)>,
}

/// Result of [`convex_hull`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ConvexHull {
    /// Elements in alphabetical order.
    pub elements: Vec<String>,
    /// Reference chemical potentials in eV/atom.
    pub references: BTreeMap<String, f64>,
    /// Analysed points in input order.
    pub points: Vec<PhasePoint>,
}

/// Computes formation energies and the distance to the lower convex hull.
///
/// When `references` is `None`, the energy per atom of each pure-element
/// entry is used as the chemical potential. The hull is solved with a
/// two-phase simplex, so any number of elements is supported.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty list, non-finite energies,
/// empty compositions or a missing elemental reference.
pub fn convex_hull(
    entries: &[PhaseEntry],
    references: Option<&BTreeMap<String, f64>>,
) -> Result<ConvexHull, Error> {
    if entries.is_empty() {
        return Err(Error::InvalidInput(
            "the convex hull needs at least one phase".to_owned(),
        ));
    }
    let mut elements = BTreeSet::new();
    for entry in entries {
        if entry.composition.is_empty() {
            return Err(Error::InvalidInput(format!(
                "phase {} has an empty composition",
                entry.name
            )));
        }
        if !entry.energy.is_finite() {
            return Err(Error::InvalidInput(format!(
                "phase {} has a non-finite energy",
                entry.name
            )));
        }
        for (element, amount) in &entry.composition {
            if !amount.is_finite() || *amount < 0.0 {
                return Err(Error::InvalidInput(format!(
                    "phase {} has an invalid amount of {element}",
                    entry.name
                )));
            }
            if *amount > 0.0 {
                elements.insert(element.clone());
            }
        }
    }
    let elements: Vec<String> = elements.into_iter().collect();
    let reference_map = match references {
        Some(map) => {
            let mut checked = BTreeMap::new();
            for element in &elements {
                let Some(value) = map.get(element) else {
                    return Err(Error::InvalidInput(format!(
                        "no reference energy for element {element}"
                    )));
                };
                if !value.is_finite() {
                    return Err(Error::InvalidInput(format!(
                        "reference energy for {element} is not finite"
                    )));
                }
                checked.insert(element.clone(), *value);
            }
            checked
        }
        None => {
            let mut derived = BTreeMap::new();
            for entry in entries {
                let total: f64 = entry.composition.values().sum();
                let non_zero: Vec<&String> = entry
                    .composition
                    .iter()
                    .filter(|(_, amount)| **amount > 0.0)
                    .map(|(element, _)| element)
                    .collect();
                if non_zero.len() == 1
                    && total > 0.0
                    && let Some(element) = non_zero.first()
                {
                    derived
                        .entry((*element).clone())
                        .or_insert(entry.energy / total);
                }
            }
            for element in &elements {
                if !derived.contains_key(element) {
                    return Err(Error::InvalidInput(format!(
                        "no pure-element phase to derive the chemical potential of {element}"
                    )));
                }
            }
            derived
        }
    };
    let mut fractions: Vec<BTreeMap<String, f64>> = Vec::with_capacity(entries.len());
    let mut formation = Vec::with_capacity(entries.len());
    for entry in entries {
        let total: f64 = entry.composition.values().sum();
        if total <= 0.0 {
            return Err(Error::InvalidInput(format!(
                "phase {} has a zero formula unit",
                entry.name
            )));
        }
        let mut fraction = BTreeMap::new();
        for element in &elements {
            fraction.insert(
                element.clone(),
                entry.composition.get(element).copied().unwrap_or(0.0) / total,
            );
        }
        let reference_energy: f64 = entry
            .composition
            .iter()
            .map(|(element, amount)| reference_map.get(element).copied().unwrap_or(0.0) * amount)
            .sum();
        formation.push((entry.energy - reference_energy) / total);
        fractions.push(fraction);
    }
    let mut points = Vec::with_capacity(entries.len());
    for index in 0..entries.len() {
        let target = &fractions[index];
        let (hull_energy, decomposition) = solve_hull(&fractions, &formation, target, &elements)?;
        let distance = formation[index] - hull_energy;
        points.push(PhasePoint {
            name: entries[index].name.clone(),
            composition: target.clone(),
            formation_energy: formation[index],
            hull_energy,
            distance,
            stable: distance <= 1e-6,
            decomposition,
        });
    }
    Ok(ConvexHull {
        elements,
        references: reference_map,
        points,
    })
}

fn solve_hull(
    fractions: &[BTreeMap<String, f64>],
    energies: &[f64],
    target: &BTreeMap<String, f64>,
    elements: &[String],
) -> Result<(f64, Vec<(usize, f64)>), Error> {
    let rows = elements.len().max(1);
    let columns = fractions.len();
    let mut matrix = vec![vec![0.0; columns]; rows];
    let mut vector = vec![0.0; rows];
    for (row, element) in elements.iter().enumerate() {
        for column in 0..columns {
            matrix[row][column] = fractions[column].get(element).copied().unwrap_or(0.0);
        }
        vector[row] = target.get(element).copied().unwrap_or(0.0);
    }
    if elements.is_empty() {
        matrix[0][0] = 1.0;
        vector[0] = 1.0;
    }
    let solution = simplex(&matrix, &vector, energies)?;
    let mut decomposition = Vec::new();
    let mut energy = 0.0;
    for (index, weight) in solution.iter().enumerate() {
        if *weight > 1e-9 {
            decomposition.push((index, *weight));
            energy += weight * energies[index];
        }
    }
    Ok((energy, decomposition))
}

/// Two-phase simplex for `min c w` subject to `A w = b`, `w >= 0`.
fn simplex(matrix: &[Vec<f64>], vector: &[f64], costs: &[f64]) -> Result<Vec<f64>, Error> {
    let rows = matrix.len();
    let columns = costs.len();
    if rows == 0 || columns == 0 {
        return Err(Error::InvalidInput(
            "the simplex needs a non-empty system".to_owned(),
        ));
    }
    let total = columns + rows;
    let mut tableau = vec![vec![0.0; total + 1]; rows + 1];
    let mut basis = vec![0usize; rows];
    for row in 0..rows {
        for column in 0..columns {
            tableau[row][column] = matrix[row][column];
        }
        tableau[row][columns + row] = 1.0;
        tableau[row][total] = vector[row];
        basis[row] = columns + row;
    }
    for column in 0..total {
        let mut sum = 0.0;
        for row in 0..rows {
            sum += tableau[row][column];
        }
        tableau[rows][column] = if column < columns { -sum } else { 1.0 - sum };
    }
    let mut sum = 0.0;
    for row in 0..rows {
        sum += vector[row];
    }
    tableau[rows][total] = -sum;
    pivot_until_done(&mut tableau, &mut basis, rows, total, Some(columns))?;
    if -tableau[rows][total] > 1e-8 {
        return Err(Error::InvalidInput(
            "the convex hull is infeasible".to_owned(),
        ));
    }
    for column in 0..total {
        let mut sum = 0.0;
        for (row, &index) in basis.iter().enumerate() {
            let cost = if index < columns { costs[index] } else { 0.0 };
            sum += cost * tableau[row][column];
        }
        tableau[rows][column] = if column < columns {
            costs[column] - sum
        } else {
            1e9
        };
    }
    let mut sum = 0.0;
    for (row, &index) in basis.iter().enumerate() {
        if index < columns {
            sum += costs[index] * tableau[row][total];
        }
    }
    tableau[rows][total] = -sum;
    pivot_until_done(&mut tableau, &mut basis, rows, total, Some(columns))?;
    let mut solution = vec![0.0; columns];
    for (row, &index) in basis.iter().enumerate() {
        if index < columns {
            solution[index] = tableau[row][total];
        }
    }
    Ok(solution)
}

fn pivot_until_done(
    tableau: &mut [Vec<f64>],
    basis: &mut [usize],
    rows: usize,
    total: usize,
    allowed: Option<usize>,
) -> Result<(), Error> {
    let limit = allowed.unwrap_or(total);
    for _ in 0..10_000 {
        let mut entering = None;
        let mut best = -1e-10;
        for column in 0..limit {
            let value = tableau[rows][column];
            if value < best {
                best = value;
                entering = Some(column);
            }
        }
        let Some(column) = entering else {
            return Ok(());
        };
        let mut leaving = None;
        let mut ratio = f64::INFINITY;
        for row in 0..rows {
            let coefficient = tableau[row][column];
            if coefficient > 1e-12 {
                let candidate = tableau[row][total] / coefficient;
                if candidate < ratio - 1e-12 {
                    ratio = candidate;
                    leaving = Some(row);
                }
            }
        }
        let Some(row) = leaving else {
            return Err(Error::InvalidInput(
                "the convex-hull linear program is unbounded".to_owned(),
            ));
        };
        basis[row] = column;
        let pivot = tableau[row][column];
        for value in 0..=total {
            tableau[row][value] /= pivot;
        }
        for other in 0..=rows {
            if other == row {
                continue;
            }
            let factor = tableau[other][column];
            if factor.abs() < 1e-18 {
                continue;
            }
            for value in 0..=total {
                tableau[other][value] -= factor * tableau[row][value];
            }
        }
    }
    Err(Error::InvalidInput(
        "the convex-hull simplex did not converge".to_owned(),
    ))
}
