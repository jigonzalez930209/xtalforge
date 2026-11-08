//! The relaxed-scan driver (Task 4.5.3).
//!
//! For each value of an axis the coordinate is moved to that value, the rest
//! of the structure is relaxed under a quadratic penalty holding the
//! coordinate, and the point records both the target and the value actually
//! reached.  Each point is seeded from the previous relaxed geometry, which is
//! what makes a scan affordable; the consequence is path dependence, and
//! `achieved` plus the per-point structure are kept so it is visible.
//!
//! **Documented subset.**  Scans operate on P1 structures, one axis and one
//! coordinate at a time.  Rigid scans, two-dimensional grids, on-disk
//! checkpointing and exact holonomic constraints are follow-up work recorded
//! in Task 4.5; the penalty and its constant are explicit here.

use serde::{Deserialize, Serialize};
use xtalforge_core::Structure;

use crate::Error;
use crate::calculator::Calculator;
use crate::math::{norm, sub};
use crate::optimize::{Minimizer, OptimizeOptions, minimize};
use crate::scan::constraints::Constraint;
use crate::scan::coordinates::ScanCoordinate;
use crate::structure_util::with_frac;

/// Options for a relaxed scan.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScanOptions {
    /// Optimizer settings used at every point; the cell is normally fixed.
    pub optimizer: OptimizeOptions,
    /// Penalty constant in eV per squared coordinate unit.
    pub penalty: f64,
    /// Minimizer used at every point.
    pub minimizer: Minimizer,
}

impl Default for ScanOptions {
    fn default() -> Self {
        Self {
            optimizer: OptimizeOptions::default(),
            penalty: 50.0,
            minimizer: Minimizer::Lbfgs,
        }
    }
}

impl ScanOptions {
    /// Checks that the options are usable.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a negative or non-finite penalty,
    /// or for bad optimizer settings.
    pub fn validate(&self) -> Result<(), Error> {
        self.optimizer.validate()?;
        if !self.penalty.is_finite() || self.penalty < 0.0 {
            return Err(Error::InvalidInput(format!(
                "scan penalty must be non-negative and finite, got {}",
                self.penalty
            )));
        }
        Ok(())
    }
}

/// One coordinate and the values it is stepped through.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScanAxis {
    /// Coordinate being scanned.
    pub coordinate: ScanCoordinate,
    /// Values, in the coordinate's units.
    pub values: Vec<f64>,
}

impl ScanAxis {
    /// Builds an axis from explicit values.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty list or a non-finite
    /// value.
    pub fn new(coordinate: ScanCoordinate, values: Vec<f64>) -> Result<Self, Error> {
        if values.is_empty() {
            return Err(Error::InvalidInput(
                "a scan axis needs at least one value".to_owned(),
            ));
        }
        if values.iter().any(|value| !value.is_finite()) {
            return Err(Error::InvalidInput("scan values must be finite".to_owned()));
        }
        Ok(Self { coordinate, values })
    }

    /// `points` evenly spaced values from `start` to `stop`, inclusive.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for fewer than two points or
    /// non-finite endpoints.
    pub fn linspace(
        coordinate: ScanCoordinate,
        start: f64,
        stop: f64,
        points: usize,
    ) -> Result<Self, Error> {
        if points < 2 {
            return Err(Error::InvalidInput(
                "a scan needs at least two points".to_owned(),
            ));
        }
        if !start.is_finite() || !stop.is_finite() {
            return Err(Error::InvalidInput(
                "scan endpoints must be finite".to_owned(),
            ));
        }
        let step = (stop - start) / (points - 1) as f64;
        let values = (0..points)
            .map(|index| start + step * index as f64)
            .collect();
        Self::new(coordinate, values)
    }

    /// Coordinate label.
    #[must_use]
    pub fn label(&self) -> String {
        self.coordinate.label()
    }
}

/// One relaxed point of a scan.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScanPoint {
    /// Requested coordinate value.
    pub target: f64,
    /// Coordinate value after relaxation.
    pub achieved: f64,
    /// Physical energy in eV, with the penalty excluded.
    pub energy: f64,
    /// Largest force on a free site after relaxation, in eV/Angstrom.
    pub max_force: f64,
    /// Optimizer iterations at this point.
    pub steps: usize,
    /// Whether the optimizer reported convergence.
    pub converged: bool,
    /// Relaxed structure.
    pub structure: Structure,
}

/// Every point of a scan.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScanResult {
    /// Coordinate that was scanned.
    pub coordinate: ScanCoordinate,
    /// Points in the order they were walked.
    pub points: Vec<ScanPoint>,
}

impl ScanResult {
    /// Energies in point order, in eV.
    #[must_use]
    pub fn energies(&self) -> Vec<f64> {
        self.points.iter().map(|point| point.energy).collect()
    }

    /// The lowest finished point.
    #[must_use]
    pub fn minimum(&self) -> Option<&ScanPoint> {
        self.points
            .iter()
            .filter(|point| point.energy.is_finite())
            .min_by(|first, second| first.energy.total_cmp(&second.energy))
    }
}

/// Relaxes the structure at every value of one axis.
///
/// `calculator` must evaluate the same structure the scan moves; the returned
/// points carry the relaxed structures for the next step and for the caller.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-P1 structure, invalid options,
/// an index outside the structure, or an invalid coordinate, and any
/// calculator error raised at a point.
pub fn relaxed_scan(
    calculator: &dyn Calculator,
    structure: &Structure,
    axis: &ScanAxis,
    options: &ScanOptions,
) -> Result<ScanResult, Error> {
    options.validate()?;
    if !structure.is_p1() {
        return Err(Error::InvalidInput(
            "relaxed scans operate on P1 structures; reduce the structure to P1 first".to_owned(),
        ));
    }
    axis.coordinate.validate(structure.sites().len())?;
    if axis.values.is_empty() {
        return Err(Error::InvalidInput(
            "a scan axis needs at least one value".to_owned(),
        ));
    }

    let mut current = structure.clone();
    let mut points = Vec::with_capacity(axis.values.len());
    for target in &axis.values {
        let lattice = *current.lattice();
        let mut positions = current.cart_sites();
        axis.coordinate.set(&mut positions, &lattice, *target)?;
        let frac: Vec<[f64; 3]> = positions
            .iter()
            .map(|position| lattice.to_frac(*position))
            .collect();
        let seeded = with_frac(&current, lattice, &frac)?;
        let constraint = Constraint::new(axis.coordinate.clone(), *target)?;
        let penalized = PenalizedCalculator {
            inner: calculator,
            constraints: vec![constraint],
            constant: options.penalty,
        };
        let result = minimize(&penalized, &seeded, &options.optimizer, options.minimizer)?;
        let achieved = axis
            .coordinate
            .value(&result.structure.cart_sites(), result.structure.lattice())?;
        let energy = calculator.energy(&result.structure)?;
        let (_, forces) = calculator.energy_forces(&result.structure)?;
        let max_force = forces
            .iter()
            .map(|force| norm(*force))
            .fold(0.0_f64, f64::max);
        points.push(ScanPoint {
            target: *target,
            achieved,
            energy,
            max_force,
            steps: result.steps,
            converged: result.converged,
            structure: result.structure.clone(),
        });
        current = result.structure;
    }
    Ok(ScanResult {
        coordinate: axis.coordinate.clone(),
        points,
    })
}

/// A calculator with a quadratic penalty added for each constraint.
struct PenalizedCalculator<'a> {
    inner: &'a dyn Calculator,
    constraints: Vec<Constraint>,
    constant: f64,
}

impl Calculator for PenalizedCalculator<'_> {
    fn name(&self) -> &'static str {
        "penalized"
    }

    fn energy(&self, structure: &Structure) -> Result<f64, Error> {
        let mut energy = self.inner.energy(structure)?;
        let positions = structure.cart_sites();
        for constraint in &self.constraints {
            let (penalty, _) =
                constraint.penalty(&positions, structure.lattice(), self.constant)?;
            energy += penalty;
        }
        Ok(energy)
    }

    fn energy_forces(&self, structure: &Structure) -> Result<(f64, Vec<[f64; 3]>), Error> {
        let (mut energy, mut forces) = self.inner.energy_forces(structure)?;
        let positions = structure.cart_sites();
        for constraint in &self.constraints {
            let (penalty, gradient) =
                constraint.penalty(&positions, structure.lattice(), self.constant)?;
            energy += penalty;
            for (force, component) in forces.iter_mut().zip(&gradient) {
                *force = sub(*force, *component);
            }
        }
        Ok((energy, forces))
    }

    fn provides_forces(&self) -> bool {
        true
    }

    fn provides_stress(&self) -> bool {
        false
    }
}
