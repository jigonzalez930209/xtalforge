//! YAML workflow reading (Task 17.2.1).
//!
//! The YAML subset is exercised through the public `Workflow` API only: the
//! same fixture as JSON must parse to the same value, and every syntax error
//! must carry line and column context.

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::path::PathBuf;

use xtalforge_agent::workflow::{StepKind, Workflow};

const JSON_FIXTURE: &str = include_str!("fixtures/workflow.json");
const YAML_FIXTURE: &str = include_str!("fixtures/workflow.yaml");

fn temp_dir(tag: &str) -> PathBuf {
    let path = std::env::temp_dir().join(format!("xtalforge-wf-yaml-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&path);
    std::fs::create_dir_all(&path).expect("temp dir");
    path
}

#[test]
fn the_yaml_fixture_matches_the_json_fixture() {
    let json = Workflow::from_json(JSON_FIXTURE).expect("JSON fixture");
    let yaml = Workflow::from_yaml(YAML_FIXTURE).expect("YAML fixture");
    assert_eq!(yaml, json);
    assert_eq!(
        yaml.validate().expect("valid YAML"),
        json.validate().expect("valid JSON")
    );
    let json_as_yaml = Workflow::from_yaml(JSON_FIXTURE).expect("JSON is valid YAML");
    assert_eq!(json_as_yaml, json);
}

#[test]
fn a_flow_style_root_mapping_parses() {
    let text = "{workflow_version: 1, name: x, steps: [{kind: read, path: '{{input}}'}]}";
    let workflow = Workflow::from_yaml(text).expect("flow root");
    assert_eq!(workflow.name, "x");
    assert_eq!(workflow.steps.len(), 1);
    assert_eq!(Workflow::from_text(text).expect("detected"), workflow);
}

#[test]
fn nested_steps_flow_collections_and_scalars_parse() {
    let text = r#"
workflow_version: 1
name: nested
description: "Nested: steps and inputs"
steps:
- kind: read
  path: '{{input}}'
- kind: transform
  name: expand
  verb: supercell
  arguments:
    repeats: [2, 2, 2]
    label: "cell: doubled"
    enabled: true
    weight: 0.5
    note: null
- kind: analyse
  name: symmetry
  verb: symmetry
  arguments: {symprec: 0.001, quiet: false}
- kind: write
  name: copy
  path: "{{dir}}/out/{{stem}}.cif"
"#;
    let workflow = Workflow::from_yaml(text).expect("nested YAML");
    let summary = workflow.validate().expect("valid");
    assert_eq!(summary.steps.len(), 4);
    assert_eq!(summary.steps[1].kind, StepKind::Transform);
    assert_eq!(summary.steps[1].verb.as_deref(), Some("supercell"));

    let transform = workflow.steps[1].arguments.as_ref().expect("arguments");
    assert_eq!(transform["repeats"], serde_json::json!([2, 2, 2]));
    assert_eq!(transform["label"], "cell: doubled");
    assert_eq!(transform["enabled"], true);
    assert_eq!(transform["weight"], 0.5);
    assert!(transform["note"].is_null());

    let analyse = workflow.steps[2].arguments.as_ref().expect("arguments");
    assert_eq!(analyse["symprec"], 0.001);
    assert_eq!(analyse["quiet"], false);
}

#[test]
fn quoted_values_keep_colons_and_comment_markers() {
    let text = r#"
# leading comment
workflow_version: 1
name: quoted # trailing comment
description: "Inspect: every input # still text"  # a real comment
steps:
  - kind: read
    path: '{{dir}}/in:put.cif'  # single quotes keep the colon
  - kind: write
    path: "{{dir}}/out/{{stem}}.cif"
"#;
    let workflow = Workflow::from_yaml(text).expect("quoted YAML");
    assert_eq!(workflow.name, "quoted");
    assert_eq!(
        workflow.description.as_deref(),
        Some("Inspect: every input # still text")
    );
    assert_eq!(
        workflow.steps[0].path.as_deref(),
        Some("{{dir}}/in:put.cif")
    );
}

#[test]
fn unquoted_template_placeholders_stay_plain_scalars() {
    let text = "workflow_version: 1\nname: x\nsteps:\n  - kind: read\n    path: {{input}}\n";
    let workflow = Workflow::from_yaml(text).expect("template path");
    assert_eq!(
        workflow.validate().expect("valid").steps[0].path.as_deref(),
        Some("{{input}}")
    );
}

#[test]
fn invalid_yaml_reports_a_line_and_column() {
    let text = "workflow_version: 1\nname: \"broken\nsteps:\n  - kind: read\n";
    let error = Workflow::from_yaml(text).expect_err("unterminated quote");
    assert_eq!(error.code(), "workflow_parse");
    let message = error.to_string();
    assert!(message.contains("line 2"), "{message}");
    assert!(message.contains("column"), "{message}");

    let text = "workflow_version: 1\nname: x\nsteps:\n  - kind: read\n   path: oops.cif\n";
    let error = Workflow::from_yaml(text).expect_err("bad indentation");
    assert_eq!(error.code(), "workflow_parse");
    assert!(error.to_string().contains("line 5"), "{error}");
}

#[test]
fn unknown_fields_still_fail_validation() {
    let text = "workflow_version: 1\nname: x\nsteps:\n  - kind: read\nsurprise: true\n";
    let error = Workflow::from_yaml(text).expect_err("unknown field");
    assert_eq!(error.code(), "workflow_parse");
    assert!(error.to_string().contains("unknown field"), "{error}");
}

#[test]
fn unsupported_yaml_features_are_rejected() {
    let cases = [
        (
            "workflow_version: 1\nname: &name anchored\nsteps:\n  - kind: read\n",
            "anchors",
        ),
        (
            "workflow_version: 1\nname: x\nsteps:\n  - kind: read\n    path: *path\n",
            "aliases",
        ),
        (
            "workflow_version: 1\nname: x\ndescription: !text hello\nsteps:\n  - kind: read\n",
            "tags",
        ),
        (
            "workflow_version: 1\nname: x\nsteps:\n  - kind: read\n---\nname: y\n",
            "multi-document",
        ),
        (
            "workflow_version: 1\nname: x\ndescription: |\n  two lines\nsteps:\n  - kind: read\n",
            "block scalars",
        ),
    ];
    for (text, needle) in cases {
        let error = Workflow::from_yaml(text).expect_err(needle);
        assert_eq!(error.code(), "workflow_parse");
        assert!(
            error.to_string().contains(needle),
            "expected {needle:?} in {error}"
        );
        assert!(error.to_string().contains("line"), "{error}");
    }
}

#[test]
fn from_text_detects_json_and_yaml() {
    let json = Workflow::from_text(JSON_FIXTURE).expect("JSON text");
    let yaml = Workflow::from_text(YAML_FIXTURE).expect("YAML text");
    assert_eq!(json, yaml);
    let marked = format!("---\n{YAML_FIXTURE}");
    assert_eq!(Workflow::from_text(&marked).expect("marker"), yaml);
    assert_eq!(Workflow::from_yaml(&marked).expect("marker"), yaml);
}

#[test]
fn from_path_selects_yaml_by_extension() {
    let directory = temp_dir("path");
    let yaml_path = directory.join("workflow.yaml");
    let yml_path = directory.join("workflow.yml");
    let upper_path = directory.join("workflow.YAML");
    std::fs::write(&yaml_path, YAML_FIXTURE).expect("write yaml");
    std::fs::copy(&yaml_path, &yml_path).expect("copy yml");
    std::fs::copy(&yaml_path, &upper_path).expect("copy upper");

    let expected = Workflow::from_json(JSON_FIXTURE).expect("JSON fixture");
    assert_eq!(Workflow::from_path(&yaml_path).expect("yaml"), expected);
    assert_eq!(Workflow::from_path(&yml_path).expect("yml"), expected);
    assert_eq!(Workflow::from_path(&upper_path).expect("YAML"), expected);

    let json_path = directory.join("workflow.json");
    std::fs::write(&json_path, JSON_FIXTURE).expect("write json");
    assert_eq!(Workflow::from_path(&json_path).expect("json"), expected);

    let plain_path = directory.join("workflow");
    std::fs::copy(&yaml_path, &plain_path).expect("copy extensionless");
    assert_eq!(Workflow::from_path(&plain_path).expect("auto"), expected);

    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn from_path_reports_io_errors() {
    let directory = temp_dir("missing");
    let path = directory.join("missing.yaml");
    let error = Workflow::from_path(&path).expect_err("missing file");
    assert_eq!(error.code(), "io");
    let _ = std::fs::remove_dir_all(&directory);
}
