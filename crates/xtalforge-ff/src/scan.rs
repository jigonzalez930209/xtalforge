//! Coordinates, constraints and relaxed scans (Tasks 4.5 and 14.2.2/14.2.3).
//!
//! [`ScanCoordinate`] names a distance, angle or torsion by atom index, a cell
//! parameter, or a collective variable; [`Constraint`] holds one near a target
//! with a quadratic penalty; and [`relaxed_scan`] steps an axis, seeds each
//! point from the previous relaxed geometry, and records both the target and
//! the value reached.  Scans are P1 and one axis at a time.
//!
//! [`relaxed_scan_with`] adds cancellation, progress, trajectory recording and
//! an on-disk [`ScanCheckpoint`] that resumes a scan at its first unfinished
//! point.  The exact holonomic projection and rigid scans remain follow-up
//! work recorded in the roadmap.

pub mod checkpoint;
pub mod constraints;
pub mod coordinates;
pub mod driver;

mod cell_parameter;
mod collective;
mod internal;
mod penalty;

pub use cell_parameter::CellParameter;
pub use checkpoint::ScanCheckpoint;
pub use collective::CollectiveVariable;
pub use constraints::Constraint;
pub use coordinates::{ScanCoordinate, ScanGeometry};
pub use driver::{ScanAxis, ScanOptions, ScanPoint, ScanResult, relaxed_scan, relaxed_scan_with};
