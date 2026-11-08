//! Integration tests for reports and the test stub module (Task 5.10).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::json;
use xtalforge_modules::{
    Module, ModuleInfo, ModuleRegistry, Report, ReportSection, TestModule, test_module,
};

#[test]
fn a_report_renders_to_markdown() {
    let report = Report {
        run_id: "run-7".to_owned(),
        title: "PXRD of Si".to_owned(),
        summary: "One phase, three reflections.".to_owned(),
        sections: vec![ReportSection {
            heading: "Peaks".to_owned(),
            body: "The three strongest peaks are 111, 220 and 311.".to_owned(),
        }],
        references: vec!["Crystal-Builder golden data".to_owned()],
    };
    let markdown = report.to_markdown();
    assert!(markdown.starts_with("# PXRD of Si\n\n"));
    assert!(markdown.contains("Run: `run-7`"));
    assert!(markdown.contains("## Peaks"));
    assert!(markdown.contains("## References"));
    assert!(markdown.contains("- Crystal-Builder golden data"));
    assert!(markdown.ends_with('\n'));
    assert!(!markdown.ends_with("\n\n"));
}

#[test]
fn empty_summary_and_references_are_omitted() {
    let report = Report {
        run_id: String::new(),
        title: "Empty".to_owned(),
        summary: "   ".to_owned(),
        sections: Vec::new(),
        references: Vec::new(),
    };
    let markdown = report.to_markdown();
    assert_eq!(markdown, "# Empty\n");
}

#[test]
fn the_test_module_reports_availability_and_a_result() {
    let mut registry = ModuleRegistry::new();
    registry.register(test_module()).unwrap();
    let module = registry.get("test").unwrap();
    assert_eq!(module.info().name, "test");
    assert!(module.available());
    let result = module.run(&json!({})).unwrap();
    assert_eq!(result["ok"], json!(true));
    assert_eq!(result["module"], json!("test"));
}

#[test]
fn the_stub_type_is_usable_directly() {
    let module = TestModule;
    assert!(module.available());
    assert_eq!(module.info().name, "test");
}

#[derive(Debug)]
struct Sentinel;

impl Module for Sentinel {
    fn info(&self) -> ModuleInfo {
        ModuleInfo {
            name: "sentinel".to_owned(),
            description: "No run implementation".to_owned(),
            version: "1".to_owned(),
            parameters: Vec::new(),
        }
    }
}

#[test]
fn a_module_without_a_run_implementation_is_typed() {
    let module = Sentinel;
    assert!(module.available());
    assert!(module.run(&json!({})).is_err());
}
