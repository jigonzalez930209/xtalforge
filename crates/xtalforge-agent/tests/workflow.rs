//! Workflow model, schema and input expansion (Tasks 16.1.1 and 16.1.2).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use std::path::PathBuf;

use serde_json::json;
use xtalforge_agent::workflow::{
    StepKind, TemplateContext, WORKFLOW_SCHEMA_ID, WORKFLOW_VERSION, Workflow, expand_inputs,
    render_template, schema_json, schema_text, tool_schema,
};

const FIXTURE: &str = include_str!("fixtures/workflow.json");

fn temp_dir(tag: &str) -> PathBuf {
    let path =
        std::env::temp_dir().join(format!("xtalforge-wf-model-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&path);
    std::fs::create_dir_all(&path).expect("temp dir");
    path
}

#[test]
fn the_fixture_parses_and_summarizes() {
    let workflow = Workflow::from_json(FIXTURE).expect("fixture");
    let summary = workflow.validate().expect("valid");
    assert_eq!(summary.workflow_version, WORKFLOW_VERSION);
    assert_eq!(summary.name, "inspect-copy");
    assert_eq!(summary.steps.len(), 3);
    assert_eq!(summary.steps[0].kind, StepKind::Read);
    assert_eq!(summary.steps[0].path.as_deref(), Some("{{input}}"));
    assert_eq!(summary.steps[1].kind, StepKind::Analyse);
    assert_eq!(summary.steps[1].verb.as_deref(), Some("info"));
    assert_eq!(summary.steps[2].kind, StepKind::Write);
    assert_eq!(
        summary.steps[2].path.as_deref(),
        Some("{{dir}}/out/{{stem}}.cif")
    );
}

#[test]
fn a_compute_step_names_a_read_only_verb() {
    let text = json!({
        "workflow_version": 1,
        "name": "energy",
        "steps": [
            {"kind": "read"},
            {"kind": "compute", "verb": "uff_energy"},
            {"kind": "transform", "verb": "supercell", "arguments": {"repeats": [2, 2, 2]}}
        ]
    });
    let workflow = Workflow::from_value(text).expect("valid document");
    let summary = workflow.validate().expect("valid");
    assert_eq!(summary.steps[1].kind, StepKind::Compute);
}

#[test]
fn malformed_json_is_a_typed_parse_error() {
    let error = Workflow::from_json("{ not json").expect_err("must fail");
    assert_eq!(error.code(), "workflow_parse");
    assert_eq!(error.category(), "input");
}

#[test]
fn an_unknown_field_is_a_parse_error() {
    let text = json!({
        "workflow_version": 1,
        "name": "x",
        "steps": [{"kind": "read"}],
        "surprise": true
    });
    let error = Workflow::from_value(text).expect_err("must fail");
    assert_eq!(error.code(), "workflow_parse");
}

#[test]
fn another_version_is_a_version_error() {
    let text = json!({
        "workflow_version": 2,
        "name": "x",
        "steps": [{"kind": "read"}]
    });
    let error = Workflow::from_value(text).expect_err("must fail");
    assert_eq!(error.code(), "workflow_version");
    assert!(error.to_string().contains("version 2"));
}

#[test]
fn the_first_step_must_read() {
    let text = json!({
        "workflow_version": 1,
        "name": "x",
        "steps": [{"kind": "analyse", "verb": "info"}]
    });
    let error = Workflow::from_value(text).expect_err("must fail");
    assert_eq!(error.code(), "workflow_invalid");
    assert_eq!(error.path(), Some("steps[0].kind"));
}

#[test]
fn mutation_and_kind_must_agree() {
    let text = json!({
        "workflow_version": 1,
        "name": "x",
        "steps": [
            {"kind": "read"},
            {"kind": "transform", "verb": "info"}
        ]
    });
    let error = Workflow::from_value(text).expect_err("must fail");
    assert_eq!(error.path(), Some("steps[1].kind"));

    let text = json!({
        "workflow_version": 1,
        "name": "x",
        "steps": [
            {"kind": "read"},
            {"kind": "analyse", "verb": "supercell", "arguments": {"repeats": [2, 2, 2]}}
        ]
    });
    let error = Workflow::from_value(text).expect_err("must fail");
    assert_eq!(error.path(), Some("steps[1].kind"));
}

#[test]
fn unknown_and_plumbing_verbs_are_rejected() {
    let text = json!({
        "workflow_version": 1,
        "name": "x",
        "steps": [
            {"kind": "read"},
            {"kind": "analyse", "verb": "teleport"}
        ]
    });
    let error = Workflow::from_value(text).expect_err("must fail");
    assert_eq!(error.path(), Some("steps[1].verb"));

    let text = json!({
        "workflow_version": 1,
        "name": "x",
        "steps": [
            {"kind": "read"},
            {"kind": "analyse", "verb": "load", "arguments": {"path": "a.cif"}}
        ]
    });
    let error = Workflow::from_value(text).expect_err("must fail");
    assert!(error.to_string().contains("session plumbing"));
}

#[test]
fn read_and_write_rules_are_enforced() {
    let text = json!({
        "workflow_version": 1,
        "name": "x",
        "steps": [{"kind": "read", "verb": "info"}]
    });
    let error = Workflow::from_value(text).expect_err("must fail");
    assert_eq!(error.path(), Some("steps[0].verb"));

    let text = json!({
        "workflow_version": 1,
        "name": "x",
        "steps": [{"kind": "read"}, {"kind": "write"}]
    });
    let error = Workflow::from_value(text).expect_err("must fail");
    assert_eq!(error.path(), Some("steps[1].path"));

    let text = json!({
        "workflow_version": 1,
        "name": "x",
        "steps": [{"kind": "write", "path": "out.cif"}]
    });
    let error = Workflow::from_value(text).expect_err("must fail");
    assert_eq!(error.path(), Some("steps[0].kind"));
}

#[test]
fn formats_and_placeholders_are_checked() {
    let text = json!({
        "workflow_version": 1,
        "name": "x",
        "steps": [{"kind": "read", "format": "nope"}]
    });
    let error = Workflow::from_value(text).expect_err("must fail");
    assert_eq!(error.path(), Some("steps[0].format"));

    let text = json!({
        "workflow_version": 1,
        "name": "x",
        "steps": [{"kind": "read"}, {"kind": "write", "path": "{{nope}}/x.cif"}]
    });
    let error = Workflow::from_value(text).expect_err("must fail");
    assert_eq!(error.path(), Some("steps[1].path"));
    assert!(error.to_string().contains("nope"));
}

#[test]
fn step_names_are_unique_and_safe() {
    let text = json!({
        "workflow_version": 1,
        "name": "x",
        "steps": [
            {"kind": "read", "name": "a"},
            {"kind": "analyse", "name": "a", "verb": "info"}
        ]
    });
    let error = Workflow::from_value(text).expect_err("must fail");
    assert!(error.to_string().contains("duplicate step name"));

    let text = json!({
        "workflow_version": 1,
        "name": "x",
        "steps": [{"kind": "read", "name": "bad name"}]
    });
    let error = Workflow::from_value(text).expect_err("must fail");
    assert_eq!(error.path(), Some("steps[0].name"));
}

#[test]
fn templates_render_every_placeholder() {
    let context = TemplateContext::new(std::path::Path::new("/tmp/run/NaCl.cif"), 3);
    let rendered =
        render_template("{{dir}}/{{index}}/{{stem}}-{{ext}}.{{ext}}", &context).expect("render");
    assert_eq!(rendered, "/tmp/run/3/NaCl-cif.cif");
    assert_eq!(
        render_template("{{input}}", &context).expect("render"),
        "/tmp/run/NaCl.cif"
    );
    let error = render_template("{{nope}}", &context).expect_err("must fail");
    assert_eq!(error.code(), "workflow_invalid");
}

#[test]
fn the_document_round_trips() {
    let workflow = Workflow::from_json(FIXTURE).expect("fixture");
    let text = workflow.to_json_pretty().expect("serialize");
    let again = Workflow::from_json(&text).expect("reparse");
    assert_eq!(workflow, again);
}

#[test]
fn the_published_schema_describes_the_model() {
    let schema = schema_json();
    assert_eq!(schema["$id"], WORKFLOW_SCHEMA_ID);
    assert_eq!(
        schema["$schema"],
        "https://json-schema.org/draft/2020-12/schema"
    );
    assert_eq!(schema["properties"]["workflow_version"]["const"], 1);
    assert_eq!(schema["$defs"]["step"]["required"][0], "kind");
    assert_eq!(schema["x-xtalforge"]["schema_version"], WORKFLOW_VERSION);
    assert_eq!(schema["x-xtalforge"]["yaml"], false);
    let presets = schema["x-xtalforge"]["presets"]
        .as_array()
        .expect("presets");
    assert!(presets.iter().any(|value| value == "pxrd"));
    let placeholders = schema["x-xtalforge"]["placeholders"]
        .as_array()
        .expect("placeholders");
    assert!(placeholders.iter().any(|value| value == "input"));
    let transform = schema["x-xtalforge"]["verbs"]["transform"]
        .as_array()
        .expect("transform verbs");
    assert!(transform.iter().any(|value| value == "supercell"));
    assert!(!transform.iter().any(|value| value == "info"));
    let compute = schema["x-xtalforge"]["verbs"]["compute"]
        .as_array()
        .expect("compute verbs");
    assert!(compute.iter().any(|value| value == "uff_energy"));
    let text = schema_text().expect("schema text");
    assert!(text.contains("workflow_version"));
}

#[test]
fn the_workflow_tools_advertise_their_arguments() {
    let run = tool_schema("workflow_run").expect("workflow_run schema");
    assert_eq!(run["required"][0], "inputs");
    assert_eq!(run["properties"]["inputs"]["items"]["type"], "string");
    assert_eq!(run["properties"]["jobs"]["type"], "integer");
    let validate = tool_schema("workflow_validate").expect("validate schema");
    assert!(validate["properties"]["workflow"].is_object());
    assert!(tool_schema("info").is_none());
}

#[test]
fn globs_expand_sorted_and_deduplicated() {
    let directory = temp_dir("glob");
    std::fs::write(directory.join("b.cif"), "b").expect("write");
    std::fs::write(directory.join("a.cif"), "a").expect("write");
    std::fs::write(directory.join("notes.txt"), "x").expect("write");
    let nested = directory.join("nested");
    std::fs::create_dir_all(&nested).expect("nested");
    std::fs::write(nested.join("c.cif"), "c").expect("write");

    let pattern = format!("{}/*.cif", directory.display());
    let top = expand_inputs([pattern.as_str()]).expect("expand");
    assert_eq!(top, vec![directory.join("a.cif"), directory.join("b.cif")]);

    let deep = format!("{}/**/*.cif", directory.display());
    let all = expand_inputs([deep.as_str(), pattern.as_str()]).expect("expand");
    assert_eq!(all.len(), 3);
    assert_eq!(all[2], nested.join("c.cif"));

    let question = format!("{}/?.cif", directory.display());
    let one = expand_inputs([question.as_str()]).expect("expand");
    assert_eq!(one.len(), 2);

    let literal = directory.join("nested");
    let direct = expand_inputs([literal.to_string_lossy().as_ref()]).expect("expand");
    assert_eq!(direct, vec![nested.join("c.cif")]);

    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn a_missing_literal_input_is_an_error() {
    let directory = temp_dir("missing");
    let missing = directory.join("nope.cif");
    let error = expand_inputs([missing.to_string_lossy().as_ref()]).expect_err("must fail");
    assert_eq!(error.code(), "inputs");
    assert!(error.to_string().contains("does not exist"));
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn a_yaml_file_holding_json_is_accepted() {
    let directory = temp_dir("yaml");
    let path = directory.join("workflow.yaml");
    std::fs::write(&path, FIXTURE).expect("write");
    let workflow = Workflow::from_path(&path).expect("JSON is a YAML subset");
    assert_eq!(workflow.name, "inspect-copy");

    std::fs::write(&path, "workflow_version: 1\nname: x\nsteps: []\n").expect("write");
    let error = Workflow::from_path(&path).expect_err("YAML-only syntax is refused");
    assert_eq!(error.code(), "workflow_parse");
    assert!(error.to_string().contains("JSON only"));
    let _ = std::fs::remove_dir_all(&directory);
}
