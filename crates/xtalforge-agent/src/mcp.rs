//! MCP server over stdio (Task 8.2.1).
//!
//! The Model Context Protocol is JSON-RPC 2.0 with newline-delimited messages
//! on stdin/stdout. This is a minimal, pure-Rust implementation: one tool per
//! verb, named as the verb, with a JSON Schema built from the shared
//! catalogue. The choice of a local implementation over a crate is recorded
//! in `docs/algorithms/agent.md`.
//!
//! Supported methods: `initialize`, `notifications/initialized`, `ping`,
//! `tools/list`, `tools/call` and `notifications/cancelled` (accepted and
//! ignored; verbs are synchronous, so a client that must interrupt one
//! terminates the process).

use std::io::{BufRead, Write};

use serde_json::{Value, json};

use crate::verbs::{ParameterInfo, ParameterKind, VerbInfo};
use crate::{Error, Session};

/// Protocol version this server implements.
pub const PROTOCOL_VERSION: &str = "2024-11-05";

/// Server name reported by `initialize`.
pub const SERVER_NAME: &str = "xtalforge";

fn schema_type(kind: ParameterKind) -> &'static str {
    match kind {
        ParameterKind::String => "string",
        ParameterKind::Integer => "integer",
        ParameterKind::Number => "number",
        ParameterKind::Boolean => "boolean",
        ParameterKind::Array => "array",
        ParameterKind::Object => "object",
    }
}

fn parameter_schema(parameter: &ParameterInfo) -> Value {
    let mut schema = json!({
        "type": schema_type(parameter.kind),
        "description": parameter.summary,
    });
    if parameter.kind == ParameterKind::Array
        && let Some(object) = schema.as_object_mut()
    {
        object.insert("items".to_owned(), json!({"type": "integer"}));
        object.insert("minItems".to_owned(), json!(3));
        object.insert("maxItems".to_owned(), json!(3));
    }
    if let Some(default) = parameter.default
        && let Ok(value) = serde_json::from_str::<Value>(default)
        && let Some(object) = schema.as_object_mut()
    {
        object.insert("default".to_owned(), value);
    }
    schema
}

/// JSON Schema for one verb's arguments.
#[must_use]
pub fn input_schema(verb: &VerbInfo) -> Value {
    let mut properties = serde_json::Map::new();
    let mut required = Vec::new();
    for parameter in verb.parameters {
        properties.insert(parameter.name.to_owned(), parameter_schema(parameter));
        if parameter.required {
            required.push(Value::String(parameter.name.to_owned()));
        }
    }
    let mut schema = json!({
        "type": "object",
        "properties": properties,
        "additionalProperties": false,
    });
    if !required.is_empty()
        && let Some(object) = schema.as_object_mut()
    {
        object.insert("required".to_owned(), Value::Array(required));
    }
    schema
}

fn tool(verb: &VerbInfo) -> Value {
    json!({
        "name": verb.name,
        "description": verb.summary,
        "inputSchema": input_schema(verb),
    })
}

fn success(id: Value, result: Value) -> Value {
    json!({"jsonrpc": "2.0", "id": id, "result": result})
}

fn failure(id: Value, code: i64, message: &str) -> Value {
    json!({"jsonrpc": "2.0", "id": id, "error": {"code": code, "message": message}})
}

/// One MCP server around one session.
#[derive(Debug)]
pub struct McpServer {
    session: Session,
}

impl Default for McpServer {
    fn default() -> Self {
        Self::new()
    }
}

impl McpServer {
    /// Creates a server with an empty session.
    #[must_use]
    pub fn new() -> Self {
        Self {
            session: Session::new(),
        }
    }

    /// The session the tools run against.
    #[must_use]
    pub fn session(&self) -> &Session {
        &self.session
    }

    /// Handles one JSON-RPC message.
    ///
    /// Returns `None` for a notification (messages without an `id`) and for
    /// methods that need no answer.
    #[must_use]
    pub fn handle(&mut self, request: &Value) -> Option<Value> {
        let id = request.get("id").cloned();
        let method = request.get("method").and_then(Value::as_str).unwrap_or("");
        if method.starts_with("notifications/") || id.is_none() {
            return None;
        }
        let id = id.unwrap_or(Value::Null);
        let result = match method {
            "initialize" => Ok(json!({
                "protocolVersion": PROTOCOL_VERSION,
                "capabilities": {"tools": {"listChanged": false}},
                "serverInfo": {
                    "name": SERVER_NAME,
                    "version": xtalforge_core::version(),
                },
                "instructions": crate::context::instructions(),
            })),
            "ping" => Ok(json!({})),
            "tools/list" => Ok(json!({
                "tools": crate::verbs::catalogue().iter().map(tool).collect::<Vec<_>>(),
            })),
            "tools/call" => {
                let params = request.get("params").cloned().unwrap_or_else(|| json!({}));
                let Some(name) = params.get("name").and_then(Value::as_str) else {
                    return Some(failure(id, -32602, "tools/call needs a tool name"));
                };
                let arguments = params
                    .get("arguments")
                    .cloned()
                    .unwrap_or_else(|| json!({}));
                match self.session.run(name, &arguments) {
                    Ok(answer) => Ok(json!({
                        "content": [{
                            "type": "text",
                            "text": serde_json::to_string(&answer).unwrap_or_default(),
                        }],
                        "isError": false,
                    })),
                    Err(error) => Ok(json!({
                        "content": [{
                            "type": "text",
                            "text": serde_json::to_string(&json!({
                                "code": error.code(),
                                "message": error.to_string(),
                            }))
                            .unwrap_or_default(),
                        }],
                        "isError": true,
                    })),
                }
            }
            other => Err((-32601_i64, format!("method not found: {other}"))),
        };
        Some(match result {
            Ok(value) => success(id, value),
            Err((code, message)) => failure(id, code, &message),
        })
    }

    /// Serves newline-delimited JSON-RPC until end of input.
    ///
    /// # Errors
    ///
    /// Returns [`Error::Io`] when reading or writing fails.
    pub fn serve<R: BufRead, W: Write>(&mut self, reader: R, mut writer: W) -> Result<(), Error> {
        for line in reader.lines() {
            let line = line.map_err(|error| Error::Io(error.to_string()))?;
            if line.trim().is_empty() {
                continue;
            }
            let request: Value = match serde_json::from_str(&line) {
                Ok(value) => value,
                Err(error) => {
                    let response = failure(Value::Null, -32700, &format!("parse error: {error}"));
                    writeln!(writer, "{response}").map_err(|error| Error::Io(error.to_string()))?;
                    continue;
                }
            };
            if let Some(response) = self.handle(&request) {
                writeln!(writer, "{response}").map_err(|error| Error::Io(error.to_string()))?;
                writer
                    .flush()
                    .map_err(|error| Error::Io(error.to_string()))?;
            }
        }
        Ok(())
    }
}

/// Serves MCP over the process standard input and output.
///
/// # Errors
///
/// Returns [`Error::Io`] when reading or writing fails.
pub fn serve_stdio() -> Result<(), Error> {
    let stdin = std::io::stdin();
    let stdout = std::io::stdout();
    McpServer::new().serve(stdin.lock(), stdout.lock())
}
