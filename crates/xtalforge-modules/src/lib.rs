//! Module registry, run records and the on-disk workspace (Phase 3, Tasks
//! 3.12 and 3.13).

pub mod builtin;
pub mod jobs;
pub mod registry;
pub mod report;
pub mod stub;
pub mod workspace;

pub use builtin::registry as builtin_registry;
pub use jobs::{Job, JobQueue, JobStatus};
pub use registry::{Module, ModuleInfo, ModuleListing, ModuleRegistry};
pub use report::{Report, ReportSection};
pub use stub::{TestModule, test_module};
pub use workspace::{RunRecord, RunStatus, Workspace, archive_records, read_archived_records};
