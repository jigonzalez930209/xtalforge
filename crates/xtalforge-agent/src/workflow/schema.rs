//! The published JSON Schema for a workflow document (Task 16.1.1).
//!
//! The schema is generated from the live verb catalogue, so the list of
//! allowed verbs cannot drift from the implementation. The runner enforces
//! the same rules with typed errors; the schema is the contract a client uses
//! to build or check a document before sending it.
//!
//! [`WORKFLOW_SCHEMA_ID`] is the stable identifier (and future documentation
//! URL) of the schema. `xtalforge run --schema` prints it and the
//! `workflow_schema` verb returns it over MCP.

use serde_json::{Value, json};

use super::model::{StepKind, WORKFLOW_VERSION};
use super::rules::RESERVED_VERBS;

/// Stable `$id` of the workflow schema.
pub const WORKFLOW_SCHEMA_ID: &str =
    "https://jigonzalez930209.github.io/xtalforge/workflow.schema.json";

/// The JSON Schema of a workflow document.
#[must_use]
pub fn schema_json() -> Value {
    json!({
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": WORKFLOW_SCHEMA_ID,
        "title": "xtalforge workflow",
        "description": "A versioned pipeline of read, transform, analyse, compute and write steps applied to every input structure.",
        "type": "object",
        "additionalProperties": false,
        "required": ["workflow_version", "name", "steps"],
        "properties": {
            "workflow_version": {
                "type": "integer",
                "const": WORKFLOW_VERSION,
                "description": "Workflow model version. This build reads version 1.",
            },
            "name": {
                "type": "string",
                "minLength": 1,
                "description": "Short name used in reports and checkpoints.",
            },
            "description": {
                "type": "string",
                "description": "Optional prose for humans and assistants.",
            },
            "steps": {
                "type": "array",
                "minItems": 1,
                "description": "Ordered pipeline. The first step must be a read.",
                "items": {"$ref": "#/$defs/step"},
            },
        },
        "$defs": {
            "step": step_schema(),
        },
        "x-xtalforge": {
            "schema_version": WORKFLOW_VERSION,
            "placeholders": ["input", "stem", "ext", "dir", "index"],
            "verbs": {
                "transform": usable_verbs(StepKind::Transform),
                "analyse": usable_verbs(StepKind::Analyse),
                "compute": usable_verbs(StepKind::Compute),
            },
            "yaml": true,
            "presets": xtalforge_core::PRESET_NAMES,
        },
    })
}

/// The schema as indented JSON text.
///
/// # Errors
///
/// Returns [`serde_json::Error`] if the schema cannot be encoded (it always
/// can; the signature keeps the function panic-free).
pub fn schema_text() -> Result<String, serde_json::Error> {
    serde_json::to_string_pretty(&schema_json())
}

fn step_schema() -> Value {
    json!({
        "type": "object",
        "additionalProperties": false,
        "required": ["kind"],
        "properties": {
            "name": {
                "type": "string",
                "pattern": "^[A-Za-z0-9_.-]+$",
                "description": "Stable label used in results; defaults to step-N (1-based).",
            },
            "kind": {
                "enum": ["read", "transform", "analyse", "compute", "write"],
                "description": "What the step does.",
            },
            "verb": {
                "type": "string",
                "description": "Agent verb for transform, analyse and compute steps.",
            },
            "arguments": {
                "type": "object",
                "description": "Verb arguments, passed through unchanged.",
            },
            "path": {
                "type": "string",
                "minLength": 1,
                "description": "Path template for read and write steps; placeholders use double braces.",
            },
            "format": {
                "type": "string",
                "description": "Registered format name or extension for read and write steps.",
            },
        },
        "allOf": [
            {
                "if": {"properties": {"kind": {"const": "read"}}, "required": ["kind"]},
                "then": {"not": {"anyOf": [{"required": ["verb"]}, {"required": ["arguments"]}]}},
            },
            {
                "if": {
                    "properties": {"kind": {"enum": ["transform", "analyse", "compute"]}},
                    "required": ["kind"],
                },
                "then": {
                    "required": ["verb"],
                    "not": {"required": ["path"]},
                },
            },
            {
                "if": {"properties": {"kind": {"const": "write"}}, "required": ["kind"]},
                "then": {
                    "required": ["path"],
                    "not": {"anyOf": [{"required": ["verb"]}, {"required": ["arguments"]}]},
                },
            },
            {
                "if": {"properties": {"kind": {"enum": ["transform", "analyse", "compute"]}}, "required": ["kind"]},
                "then": {"not": {"required": ["format"]}},
            },
        ],
    })
}

fn usable_verbs(kind: StepKind) -> Vec<&'static str> {
    crate::verbs::catalogue()
        .iter()
        .filter(|verb| !RESERVED_VERBS.contains(&verb.name))
        .filter(|verb| match kind {
            StepKind::Transform => verb.mutation,
            StepKind::Analyse | StepKind::Compute => !verb.mutation,
            StepKind::Read | StepKind::Write => false,
        })
        .map(|verb| verb.name)
        .collect()
}

/// Input schema of one of the workflow MCP tools, when the name is one.
///
/// The MCP server uses this so `workflow_*` tools advertise their real
/// argument shapes instead of the generic scalar `array` schema.
#[must_use]
pub fn tool_schema(name: &str) -> Option<Value> {
    let source = json!({
        "type": "object",
        "description": "The workflow document itself.",
        "additionalProperties": true,
    });
    let properties = match name {
        "workflow_schema" => {
            json!({
                "type": "object",
                "properties": {},
                "additionalProperties": false,
            })
        }
        "workflow_validate" => json!({
            "type": "object",
            "properties": {
                "workflow": source.clone(),
                "text": {"type": "string", "description": "Workflow JSON or YAML as text."},
                "path": {"type": "string", "description": "Workflow file to read."},
            },
            "additionalProperties": false,
        }),
        "workflow_run" => json!({
            "type": "object",
            "properties": {
                "workflow": source,
                "text": {"type": "string", "description": "Workflow JSON or YAML as text."},
                "path": {"type": "string", "description": "Workflow file to read."},
                "inputs": {
                    "type": "array",
                    "items": {"type": "string"},
                    "minItems": 1,
                    "description": "Input files, directories or glob patterns.",
                },
                "jobs": {
                    "type": "integer",
                    "minimum": 1,
                    "default": 1,
                    "description": "Items processed at once.",
                },
                "checkpoint": {
                    "type": "string",
                    "description": "Checkpoint file; an existing one is resumed.",
                },
                "resume": {
                    "type": "boolean",
                    "default": true,
                    "description": "Skip items a checkpoint already completed.",
                },
            },
            "required": ["inputs"],
            "additionalProperties": false,
        }),
        _ => return None,
    };
    Some(properties)
}
