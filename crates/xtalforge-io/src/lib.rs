//! File-format readers and writers for xtalforge.
//!
//! Phase 1 provides CIF, POSCAR and XYZ; Phase 3 adds the remaining formats,
//! the format registry, project files and trajectories. Every reader and
//! writer reports failures with a stable error code and never panics on user
//! input.
