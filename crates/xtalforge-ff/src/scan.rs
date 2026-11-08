//! Coordinates, constraints and relaxed scans (Task 4.5).
//!
//! [`ScanCoordinate`] names a distance, angle or torsion by atom index;
//! [`Constraint`] holds one near a target with a quadratic penalty; and
//! [`relaxed_scan`] steps an axis, seeds each point from the previous relaxed
//! geometry, and records both the target and the value reached.  The core
//! subset is P1, one axis and one coordinate at a time; the exact holonomic
//! projection, rigid scans and on-disk checkpointing are follow-up work
//! recorded in the roadmap.

pub mod constraints;
pub mod coordinates;
pub mod driver;

pub use constraints::Constraint;
pub use coordinates::ScanCoordinate;
pub use driver::{ScanAxis, ScanOptions, ScanPoint, ScanResult, relaxed_scan};
