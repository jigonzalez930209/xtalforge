//! Module registry, run records and the on-disk workspace (Phase 3, Tasks
//! 3.12 and 3.13).

pub mod registry;
pub mod workspace;

pub use registry::{Module, ModuleInfo, ModuleRegistry};
pub use workspace::{RunRecord, RunStatus, Workspace};
