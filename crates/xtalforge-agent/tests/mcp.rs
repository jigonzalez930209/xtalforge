//! MCP JSON-RPC server (Task 8.2.1).

#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)]

use serde_json::{Value, json};
use xtalforge_agent::mcp::{McpServer, PROTOCOL_VERSION};

const NACL: &str = include_str!("../../../tests/corpus/simple/NaCl.cif");

fn call(server: &mut McpServer, id: u64, method: &str, params: Value) -> Value {
    let request = json!({"jsonrpc": "2.0", "id": id, "method": method, "params": params});
    server.handle(&request).expect("a response")
}

#[test]
fn initialize_reports_the_protocol_and_instructions() {
    let mut server = McpServer::new();
    let response = call(&mut server, 1, "initialize", json!({}));
    assert_eq!(response["result"]["protocolVersion"], PROTOCOL_VERSION);
    assert_eq!(response["result"]["serverInfo"]["name"], "xtalforge");
    let instructions = response["result"]["instructions"].as_str().expect("text");
    assert!(instructions.contains("capabilities"));
    assert!(instructions.contains("diagnostics"));
}

#[test]
fn tools_list_matches_the_catalogue() {
    let mut server = McpServer::new();
    let response = call(&mut server, 1, "tools/list", json!({}));
    let tools = response["result"]["tools"].as_array().expect("tools");
    assert_eq!(tools.len(), xtalforge_agent::verbs::VERBS.len());
    let names: Vec<&str> = tools
        .iter()
        .filter_map(|tool| tool["name"].as_str())
        .collect();
    assert!(names.contains(&"load"));
    assert!(names.contains(&"porosity"));
    let supercell = tools
        .iter()
        .find(|tool| tool["name"] == "supercell")
        .expect("supercell tool");
    assert_eq!(supercell["inputSchema"]["type"], "object");
    assert_eq!(supercell["inputSchema"]["required"][0], "repeats");
    assert_eq!(
        supercell["inputSchema"]["properties"]["repeats"]["minItems"],
        3
    );
}

#[test]
fn a_tool_call_runs_a_verb_and_returns_its_answer() {
    let mut server = McpServer::new();
    let response = call(
        &mut server,
        1,
        "tools/call",
        json!({
            "name": "load",
            "arguments": {"text": NACL, "format": "cif", "name": "nacl"}
        }),
    );
    assert_eq!(response["result"]["isError"], false);
    let text = response["result"]["content"][0]["text"]
        .as_str()
        .expect("text");
    let answer: Value = serde_json::from_str(text).expect("answer");
    assert_eq!(answer["kind"], "loaded");
    assert_eq!(answer["document"]["formula"], "Na Cl");

    let response = call(
        &mut server,
        2,
        "tools/call",
        json!({"name": "info", "arguments": {}}),
    );
    let text = response["result"]["content"][0]["text"]
        .as_str()
        .expect("text");
    let answer: Value = serde_json::from_str(text).expect("answer");
    assert_eq!(answer["kind"], "info");
    assert!(answer["diagnostics"]["diagnostics"].is_array());
}

#[test]
fn a_failing_tool_is_an_error_result_with_the_code() {
    let mut server = McpServer::new();
    let response = call(
        &mut server,
        1,
        "tools/call",
        json!({"name": "info", "arguments": {}}),
    );
    assert_eq!(response["result"]["isError"], true);
    let text = response["result"]["content"][0]["text"]
        .as_str()
        .expect("text");
    let error: Value = serde_json::from_str(text).expect("error");
    assert_eq!(error["code"], "no_document");
}

#[test]
fn unknown_methods_and_notifications_behave() {
    let mut server = McpServer::new();
    let response = call(&mut server, 7, "resources/list", json!({}));
    assert_eq!(response["error"]["code"], -32601);
    let notification = json!({"jsonrpc": "2.0", "method": "notifications/initialized"});
    assert!(server.handle(&notification).is_none());
    let cancelled = json!({"jsonrpc": "2.0", "method": "notifications/cancelled", "params": {}});
    assert!(server.handle(&cancelled).is_none());
    let ping = call(&mut server, 8, "ping", json!({}));
    assert_eq!(ping["result"], json!({}));
}

#[test]
fn serve_reads_lines_and_reports_parse_errors() {
    let input = concat!(
        "{\"jsonrpc\":\"2.0\",\"id\":1,\"method\":\"ping\"}\n",
        "not json\n",
        "{\"jsonrpc\":\"2.0\",\"id\":2,\"method\":\"initialize\",\"params\":{}}\n"
    );
    let mut output = Vec::new();
    McpServer::new()
        .serve(std::io::Cursor::new(input.as_bytes()), &mut output)
        .expect("serve");
    let text = String::from_utf8(output).expect("utf8");
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), 3);
    let first: Value = serde_json::from_str(lines[0]).expect("first");
    assert_eq!(first["id"], 1);
    let second: Value = serde_json::from_str(lines[1]).expect("second");
    assert_eq!(second["error"]["code"], -32700);
    let third: Value = serde_json::from_str(lines[2]).expect("third");
    assert_eq!(third["id"], 2);
}
