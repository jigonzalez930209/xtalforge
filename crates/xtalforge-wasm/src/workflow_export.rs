//! In-memory workflow exports (Task 17.2.5).

use std::collections::BTreeMap;

use wasm_bindgen::prelude::*;

use crate::common::set;

/// Runs a workflow over in-memory input documents.
///
/// `workflow_text` accepts JSON or YAML; `files_json` is an object mapping
/// input paths to their text. Returns `{ report, outputs }`, where `outputs`
/// maps every written path to its text. Nothing touches a filesystem.
///
/// # Errors
///
/// Returns a JavaScript `Error` for malformed workflow/file JSON, a validation
/// failure or a failed step.
#[wasm_bindgen]
pub fn run_workflow_memory(workflow_text: &str, files_json: &str) -> Result<JsValue, JsValue> {
    let workflow = xtalforge::agent::workflow::Workflow::from_text(workflow_text)
        .map_err(|error| JsValue::from_str(&error.to_string()))?;
    let files: BTreeMap<String, String> = serde_json::from_str(files_json)
        .map_err(|error| JsValue::from_str(&format!("input map: {error}")))?;
    let result = xtalforge::agent::workflow::run_workflow_memory(&workflow, &files)
        .map_err(|error| JsValue::from_str(&error.to_string()))?;
    let outputs = js_sys::Object::new();
    for (path, text) in &result.outputs {
        set(&outputs, path, JsValue::from_str(text));
    }
    let report_json = serde_json::to_string(&result.report)
        .map_err(|error| JsValue::from_str(&format!("report: {error}")))?;
    let report = js_sys::JSON::parse(&report_json)
        .map_err(|_| JsValue::from_str("report is not valid JSON"))?;
    let object = js_sys::Object::new();
    set(&object, "report", report);
    set(&object, "outputs", outputs.into());
    Ok(object.into())
}
