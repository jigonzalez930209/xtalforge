//! Workflow tools over MCP and the capability listing (Task 16.1.3).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::{Value, json};
use xtalforge_agent::mcp::McpServer;
use xtalforge_agent::{Answer, Session};

const NACL: &str = include_str!("../../../tests/corpus/simple/NaCl.cif");
const SI: &str = include_str!("../../../tests/corpus/simple/Si.cif");
const FIXTURE: &str = include_str!("fixtures/workflow.json");

fn call(server: &mut McpServer, id: u64, method: &str, params: Value) -> Value {
    let request = json!({"jsonrpc": "2.0", "id": id, "method": method, "params": params});
    server.handle(&request).expect("a response")
}

fn text(answer: &Value) -> Value {
    serde_json::from_str(
        answer["result"]["content"][0]["text"]
            .as_str()
            .expect("text"),
    )
    .expect("answer JSON")
}

fn temp_dir(tag: &str) -> std::path::PathBuf {
    let path = std::env::temp_dir().join(format!("xtalforge-wf-mcp-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&path);
    std::fs::create_dir_all(&path).expect("temp dir");
    path
}

#[test]
fn the_catalogue_lists_the_workflow_tools() {
    let mut server = McpServer::new();
    let response = call(&mut server, 1, "tools/list", json!({}));
    let tools = response["result"]["tools"].as_array().expect("tools");
    let names: Vec<&str> = tools
        .iter()
        .filter_map(|tool| tool["name"].as_str())
        .collect();
    for name in ["workflow_schema", "workflow_validate", "workflow_run"] {
        assert!(names.contains(&name), "missing {name}");
    }
    let run = tools
        .iter()
        .find(|tool| tool["name"] == "workflow_run")
        .expect("workflow_run");
    assert_eq!(run["inputSchema"]["required"][0], "inputs");
    assert_eq!(
        run["inputSchema"]["properties"]["inputs"]["items"]["type"],
        "string"
    );
    assert_eq!(run["inputSchema"]["properties"]["jobs"]["minimum"], 1);
}

#[test]
fn capabilities_report_the_workflow_version() {
    let mut session = Session::new();
    let answer = session
        .run("capabilities", &json!({}))
        .expect("capabilities");
    let Answer::Capabilities { capabilities } = answer else {
        panic!("expected capabilities");
    };
    assert_eq!(
        capabilities.workflow_version,
        xtalforge_agent::workflow::WORKFLOW_VERSION
    );
    assert!(capabilities.verbs.iter().any(|verb| verb == "workflow_run"));
}

#[test]
fn workflow_schema_returns_the_published_schema() {
    let mut session = Session::new();
    let answer = session
        .run("workflow_schema", &json!({}))
        .expect("workflow_schema");
    let Answer::WorkflowSchema {
        workflow_version,
        schema,
    } = answer
    else {
        panic!("expected workflow_schema");
    };
    assert_eq!(workflow_version, 1);
    assert_eq!(schema["$id"], xtalforge_agent::workflow::WORKFLOW_SCHEMA_ID);
}

#[test]
fn workflow_validate_reports_valid_and_invalid_documents() {
    let mut session = Session::new();
    let answer = session
        .run("workflow_validate", &json!({"text": FIXTURE}))
        .expect("validate");
    let Answer::WorkflowValidated {
        valid,
        summary,
        problem,
    } = answer
    else {
        panic!("expected workflow_validation");
    };
    assert!(valid);
    assert!(problem.is_none());
    assert_eq!(summary.expect("summary").steps.len(), 3);

    let answer = session
        .run(
            "workflow_validate",
            &json!({"text": "{\"workflow_version\": 9}"}),
        )
        .expect("validate");
    let Answer::WorkflowValidated {
        valid,
        summary,
        problem,
    } = answer
    else {
        panic!("expected workflow_validation");
    };
    assert!(!valid);
    assert!(summary.is_none());
    assert_eq!(problem.expect("problem").code, "workflow_parse");
}

#[test]
fn workflow_run_over_mcp_runs_every_input() {
    let directory = temp_dir("run");
    std::fs::write(directory.join("NaCl.cif"), NACL).expect("write");
    std::fs::write(directory.join("Si.cif"), SI).expect("write");
    let pattern = format!("{}/*.cif", directory.display());
    let mut server = McpServer::new();
    let response = call(
        &mut server,
        1,
        "tools/call",
        json!({
            "name": "workflow_run",
            "arguments": {
                "text": FIXTURE,
                "inputs": [pattern],
                "jobs": 2,
            }
        }),
    );
    assert_eq!(response["result"]["isError"], false);
    let answer = text(&response);
    assert_eq!(answer["kind"], "workflow_run");
    assert_eq!(answer["report"]["succeeded"], 2);
    assert_eq!(answer["report"]["failed"], 0);
    assert_eq!(
        answer["report"]["items"].as_array().expect("items").len(),
        2
    );
    assert!(directory.join("out/NaCl.cif").is_file());
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn workflow_run_errors_are_isolated_and_typed() {
    let directory = temp_dir("failure");
    std::fs::write(directory.join("NaCl.cif"), NACL).expect("write");
    std::fs::write(directory.join("broken.cif"), "not a structure").expect("write");
    let pattern = format!("{}/*.cif", directory.display());
    let mut server = McpServer::new();
    let response = call(
        &mut server,
        1,
        "tools/call",
        json!({
            "name": "workflow_run",
            "arguments": {"text": FIXTURE, "inputs": [pattern]}
        }),
    );
    let answer = text(&response);
    assert_eq!(answer["report"]["failed"], 1);
    assert_eq!(answer["report"]["succeeded"], 1);

    let response = call(
        &mut server,
        2,
        "tools/call",
        json!({"name": "workflow_run", "arguments": {"text": FIXTURE}}),
    );
    assert_eq!(response["result"]["isError"], true);
    let error = text(&response);
    assert_eq!(error["code"], "invalid_argument");

    let response = call(
        &mut server,
        3,
        "tools/call",
        json!({
            "name": "workflow_run",
            "arguments": {"text": "{}", "inputs": ["*.cif"]}
        }),
    );
    assert_eq!(response["result"]["isError"], true);
    let error = text(&response);
    assert_eq!(error["code"], "workflow_parse");
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn a_workflow_file_can_be_named_instead_of_inlined() {
    let directory = temp_dir("file");
    std::fs::write(directory.join("wf.json"), FIXTURE).expect("write");
    std::fs::write(directory.join("NaCl.cif"), NACL).expect("write");
    let path = directory.join("wf.json");
    let pattern = format!("{}/*.cif", directory.display());
    let mut server = McpServer::new();
    let response = call(
        &mut server,
        1,
        "tools/call",
        json!({
            "name": "workflow_run",
            "arguments": {
                "path": path.to_string_lossy(),
                "inputs": [pattern],
                "checkpoint": directory.join("ck.json").to_string_lossy(),
            }
        }),
    );
    assert_eq!(response["result"]["isError"], false);
    let answer = text(&response);
    assert_eq!(answer["report"]["inputs"], 1);
    assert!(directory.join("ck.json").is_file());
    let _ = std::fs::remove_dir_all(&directory);
}
