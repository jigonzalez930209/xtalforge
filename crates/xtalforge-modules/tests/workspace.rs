//! Integration tests for the workspace and module registry (Tasks 3.12, 3.13).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use xtalforge_modules::{Module, ModuleInfo, ModuleRegistry, RunRecord, RunStatus, Workspace};

#[derive(Debug)]
struct Demo;

impl Module for Demo {
    fn info(&self) -> ModuleInfo {
        ModuleInfo {
            name: "demo".to_owned(),
            description: "Demo module".to_owned(),
            version: "1".to_owned(),
            parameters: vec!["x".to_owned()],
        }
    }
}

#[test]
fn registry_rejects_duplicates_and_lists_modules() {
    let mut registry = ModuleRegistry::new();
    registry.register(Box::new(Demo)).unwrap();
    assert!(registry.register(Box::new(Demo)).is_err());
    assert_eq!(registry.list().len(), 1);
    assert!(registry.get("demo").is_some());
    assert!(registry.unregister("demo"));
    assert!(registry.get("demo").is_none());
}

#[test]
fn workspace_writes_and_reads_runs() {
    let root = std::env::temp_dir().join(format!("xtalforge-workspace-{}", std::process::id()));
    let workspace = Workspace::open(&root).unwrap();
    let record = RunRecord {
        id: "run-1".to_owned(),
        module: "demo".to_owned(),
        input: "structure.cif".to_owned(),
        status: RunStatus::Done,
        outputs: vec!["out.cif".to_owned()],
        data: serde_json::json!({"energy": -1.0}),
        created: 1,
    };
    let path = workspace.write_run(&record).unwrap();
    assert!(path.exists());
    assert_eq!(workspace.read_run("run-1").unwrap(), record);
    assert_eq!(workspace.runs().unwrap(), vec!["run-1".to_owned()]);
    std::fs::remove_dir_all(&root).unwrap();
}

#[test]
fn in_memory_archive_round_trips() {
    let records = vec![RunRecord {
        id: "run-2".to_owned(),
        module: "pxrd".to_owned(),
        input: "Si.cif".to_owned(),
        status: RunStatus::Done,
        outputs: vec![],
        data: serde_json::json!({"peaks": 3}),
        created: 2,
    }];
    let bytes = xtalforge_modules::archive_records(&records).unwrap();
    assert_eq!(&bytes[..2], b"PK");
    assert_eq!(
        xtalforge_modules::read_archived_records(&bytes).unwrap(),
        records
    );
}
