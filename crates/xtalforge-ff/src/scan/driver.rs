//! The relaxed-scan driver (Task 4.5.3, extended by Tasks 14.2.2 and 14.2.3).
//!
//! For each value of an axis the coordinate is moved to that value, the rest
//! of the structure is relaxed under a quadratic penalty holding the
//! coordinate, and the point records both the target and the value actually
//! reached.  Each point is seeded from the previous relaxed geometry, which is
//! what makes a scan affordable; the consequence is path dependence, and
//! `achieved` plus the per-point structure are kept so it is visible.
//!
//! **Cell parameters.**  A [`ScanCoordinate::CellParameter`] axis imposes the
//! parameter on the seeded lattice and relaxes the sites in that cell.  When
//! the optimizer also relaxes the cell, the scanned parameter is held through
//! [`CellFreedom`](crate::CellFreedom), so the reported point is exactly at the
//! target; the other held quantities are measured from the previous point and
//! a held volume is met by the remaining free lengths.  No penalty is used for
//! a cell axis; it would have a zero Cartesian gradient.
//!
//! **Cancellation, progress and resume.**  [`relaxed_scan_with`] honours a
//! [`RunControl`], records a trajectory when
//! [`ScanOptions::trajectory`] is enabled, and writes
//! [`ScanOptions::checkpoint`] after every finished point.  A cancelled point
//! is not recorded: the checkpoint keeps the last fully relaxed geometry and a
//! re-run resumes there.  [`relaxed_scan`] is the same run without a control.
//!
//! **Documented subset.**  One axis and one coordinate at a time;
//! two-dimensional grids and exact holonomic constraints remain follow-up work.

use std::path::PathBuf;

use serde::{Deserialize, Serialize};
use xtalforge_core::Structure;

use crate::Error;
use crate::calculator::Calculator;
use crate::math::norm;
use crate::optimize::{Minimizer, OptimizeOptions, append_frames, minimize_with};
use crate::run::{ProgressEvent, RunControl};
use crate::scan::CellParameter;
use crate::scan::checkpoint::ScanCheckpoint;
use crate::scan::constraints::Constraint;
use crate::scan::coordinates::{ScanCoordinate, ScanGeometry};
use crate::scan::penalty::PenalizedCalculator;
use crate::structure_util::with_frac;
use crate::trajectory::{Trajectory, TrajectoryOptions};

/// Options for a relaxed scan.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScanOptions {
    /// Optimizer settings used at every point; the cell is normally fixed.
    pub optimizer: OptimizeOptions,
    /// Penalty constant in eV per squared coordinate unit.
    pub penalty: f64,
    /// Minimizer used at every point.
    pub minimizer: Minimizer,
    /// Per-point optimizer trajectory recording; off by default.
    ///
    /// The cadence is applied to the optimizer at every point and the frames
    /// are collected in [`ScanResult::trajectory`], tagged with the point
    /// index.  A value of `every = 0` records nothing.
    #[serde(default)]
    pub trajectory: TrajectoryOptions,
    /// Checkpoint path written after every point and read on start.
    ///
    /// A `.jsonl` path stores one line per point; anything else is one JSON
    /// document.  `None` disables checkpointing.
    #[serde(default)]
    pub checkpoint: Option<PathBuf>,
}

impl Default for ScanOptions {
    fn default() -> Self {
        Self {
            optimizer: OptimizeOptions::default(),
            penalty: 50.0,
            minimizer: Minimizer::Lbfgs,
            trajectory: TrajectoryOptions::default(),
            checkpoint: None,
        }
    }
}

impl ScanOptions {
    /// Checks that the options are usable.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for a negative or non-finite penalty,
    /// bad optimizer settings, or a checkpoint path without a file name.
    pub fn validate(&self) -> Result<(), Error> {
        self.optimizer.validate()?;
        if !self.penalty.is_finite() || self.penalty < 0.0 {
            return Err(Error::InvalidInput(format!(
                "scan penalty must be non-negative and finite, got {}",
                self.penalty
            )));
        }
        if let Some(path) = &self.checkpoint
            && path.file_name().is_none()
        {
            return Err(Error::InvalidInput(format!(
                "checkpoint path {} has no file name",
                path.display()
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

impl ScanPoint {
    /// Absolute difference between the target and the value reached.
    #[must_use]
    pub fn deviation(&self) -> f64 {
        (self.achieved - self.target).abs()
    }
}

/// Every point of a scan.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScanResult {
    /// Coordinate that was scanned.
    pub coordinate: ScanCoordinate,
    /// Points in the order they were walked.
    pub points: Vec<ScanPoint>,
    /// Whether a [`RunControl`] stopped the scan.
    #[serde(default)]
    pub cancelled: bool,
    /// Recorded optimizer frames, empty unless a cadence was set.
    #[serde(default)]
    pub trajectory: Trajectory,
    /// Points reused from a checkpoint at the start of this call.
    #[serde(default)]
    pub resumed_from: usize,
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
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-P1 structure, invalid options,
/// an index outside the structure, or a checkpoint for a different scan, and
/// any calculator error raised at a point.
pub fn relaxed_scan(
    calculator: &dyn Calculator,
    structure: &Structure,
    axis: &ScanAxis,
    options: &ScanOptions,
) -> Result<ScanResult, Error> {
    relaxed_scan_with(calculator, structure, axis, options, &mut RunControl::new())
}

/// Relaxes the structure at every value of one axis, reporting progress and
/// honouring cancellation.
///
/// # Errors
///
/// Returns the errors of [`relaxed_scan`].  Cancellation is not an error: the
/// returned result carries the finished points and `cancelled = true`.
pub fn relaxed_scan_with(
    calculator: &dyn Calculator,
    structure: &Structure,
    axis: &ScanAxis,
    options: &ScanOptions,
    control: &mut RunControl<'_>,
) -> Result<ScanResult, Error> {
    options.validate()?;
    if !structure.is_p1() {
        return Err(Error::InvalidInput(
            "relaxed scans operate on P1 structures; reduce the structure to P1 first".to_owned(),
        ));
    }
    let n_sites = structure.sites().len();
    axis.coordinate.validate(n_sites)?;
    if axis.values.is_empty() {
        return Err(Error::InvalidInput(
            "a scan axis needs at least one value".to_owned(),
        ));
    }

    let mut current = structure.clone();
    let mut points: Vec<ScanPoint> = Vec::new();
    let mut resumed_from = 0_usize;
    if let Some(path) = options.checkpoint.as_deref()
        && path.exists()
    {
        let checkpoint = ScanCheckpoint::read(path)?;
        checkpoint.matches(&axis.coordinate, &axis.values, n_sites)?;
        if let Some(last) = checkpoint.points.last() {
            current = last.structure.clone();
        }
        resumed_from = checkpoint.resumed_len();
        points = checkpoint.points;
    }

    let mut trajectory = Trajectory::new();
    for (index, target) in axis.values.iter().enumerate() {
        if index < points.len() {
            continue;
        }
        if control.is_cancelled() {
            return Ok(ScanResult {
                coordinate: axis.coordinate.clone(),
                points,
                cancelled: true,
                trajectory,
                resumed_from,
            });
        }

        let mut geometry = ScanGeometry::from_structure(&current);
        axis.coordinate.set_in(&mut geometry, *target)?;
        let mut seeded = geometry.into_structure(&current)?;

        let mut optimize = options.optimizer.clone();
        // The scan owns the cadence: frames of every point go to the scan
        // trajectory and are tagged with the point index.
        optimize.trajectory = options.trajectory;
        if let Some(parameter) = axis.coordinate.parameter() {
            optimize.cell_freedom = parameter.hold(optimize.cell_freedom);
            if optimize.cell {
                // Re-impose the cell constraints on the seed, measured from
                // the previous point: the scanned parameter keeps its new
                // target, every other held quantity comes from the previous
                // point and a held volume is met by the other free lengths.
                let mut seed_freedom = options.optimizer.cell_freedom;
                let keep_length = match parameter {
                    CellParameter::A => Some(0),
                    CellParameter::B => Some(1),
                    CellParameter::C => Some(2),
                    CellParameter::Alpha => {
                        seed_freedom = seed_freedom.with_angle(0, true);
                        None
                    }
                    CellParameter::Beta => {
                        seed_freedom = seed_freedom.with_angle(1, true);
                        None
                    }
                    CellParameter::Gamma => {
                        seed_freedom = seed_freedom.with_angle(2, true);
                        None
                    }
                };
                if let Some(index) = keep_length {
                    seed_freedom = seed_freedom.with_length(index, true);
                }
                if !seed_freedom.is_full() {
                    let lattice = seed_freedom.apply_keeping(
                        current.lattice(),
                        *seeded.lattice(),
                        keep_length,
                    )?;
                    let frac: Vec<[f64; 3]> = seeded
                        .sites()
                        .iter()
                        .map(xtalforge_core::Site::frac)
                        .collect();
                    seeded = with_frac(&seeded, lattice, &frac)?;
                }
            }
        }
        let result = if axis.coordinate.is_cell_parameter() {
            minimize_with(calculator, &seeded, &optimize, options.minimizer, control)?
        } else {
            let constraint = Constraint::new(axis.coordinate.clone(), *target)?;
            let penalized = PenalizedCalculator::new(calculator, vec![constraint], options.penalty);
            minimize_with(&penalized, &seeded, &optimize, options.minimizer, control)?
        };
        if result.cancelled {
            write_checkpoint(options, axis, n_sites, &points)?;
            return Ok(ScanResult {
                coordinate: axis.coordinate.clone(),
                points,
                cancelled: true,
                trajectory,
                resumed_from,
            });
        }

        let achieved = axis
            .coordinate
            .value(&result.structure.cart_sites(), result.structure.lattice())?;
        let energy = calculator.energy(&result.structure)?;
        let (_, forces) = calculator.energy_forces(&result.structure)?;
        let max_force = forces
            .iter()
            .map(|force| norm(*force))
            .fold(0.0_f64, f64::max);
        let steps = result.steps;
        let converged = result.converged;
        let relaxed = result.structure;
        points.push(ScanPoint {
            target: *target,
            achieved,
            energy,
            max_force,
            steps,
            converged,
            structure: relaxed.clone(),
        });
        current = relaxed;
        append_frames(&mut trajectory, result.trajectory, index);
        control.report(ProgressEvent::Scan {
            point: index,
            target: *target,
            achieved,
            energy,
        });
        write_checkpoint(options, axis, n_sites, &points)?;
    }

    Ok(ScanResult {
        coordinate: axis.coordinate.clone(),
        points,
        cancelled: false,
        trajectory,
        resumed_from,
    })
}

/// Writes the checkpoint when one is configured.
fn write_checkpoint(
    options: &ScanOptions,
    axis: &ScanAxis,
    n_sites: usize,
    points: &[ScanPoint],
) -> Result<(), Error> {
    let Some(path) = &options.checkpoint else {
        return Ok(());
    };
    let checkpoint = ScanCheckpoint {
        coordinate: axis.coordinate.clone(),
        values: axis.values.clone(),
        n_sites,
        points: points.to_vec(),
    };
    checkpoint.write(path)
}
