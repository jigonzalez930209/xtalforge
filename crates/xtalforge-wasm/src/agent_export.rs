//! WebAssembly bindings for the agent layer (Task 8.1).
//!
//! Filesystem verbs (`load` from a path, `save`, `render` to a path) are not
//! available in the browser; `load` from text and every in-memory verb work.

use wasm_bindgen::prelude::*;

use xtalforge::agent::Session;

/// Runtime capability report as JSON.
#[wasm_bindgen]
pub fn agent_capabilities() -> String {
    let report = xtalforge::agent::capability_report();
    serde_json::to_string(&report).unwrap_or_else(|_| "{}".to_owned())
}

/// The verb catalogue as JSON.
#[wasm_bindgen]
pub fn agent_verbs() -> String {
    serde_json::to_string(xtalforge::agent::catalogue()).unwrap_or_else(|_| "[]".to_owned())
}

/// Runs a JSON script against a fresh session and returns the answers.
///
/// # Errors
///
/// Returns a JavaScript object with `code` and `message` properties.
#[wasm_bindgen]
pub fn agent_run(script_json: &str) -> Result<String, JsValue> {
    let answers = xtalforge::agent::run_verbs(script_json).map_err(agent_error)?;
    serde_json::to_string(&answers).map_err(|error| JsValue::from_str(&error.to_string()))
}

fn agent_error(error: xtalforge::agent::Error) -> JsValue {
    let object = js_sys::Object::new();
    let _ = js_sys::Reflect::set(
        &object,
        &JsValue::from_str("code"),
        &JsValue::from_str(error.code()),
    );
    let _ = js_sys::Reflect::set(
        &object,
        &JsValue::from_str("message"),
        &JsValue::from_str(&error.to_string()),
    );
    object.into()
}

/// A stateful session an assistant drives in the browser.
#[wasm_bindgen]
pub struct AgentSession {
    inner: Session,
}

#[wasm_bindgen]
impl AgentSession {
    /// Creates an empty session.
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self {
            inner: Session::new(),
        }
    }

    /// Runs one verb with a JSON object and returns the answer as JSON.
    ///
    /// # Errors
    ///
    /// Returns a JavaScript object with `code` and `message` properties.
    pub fn run(&mut self, verb: &str, arguments_json: &str) -> Result<String, JsValue> {
        let arguments: serde_json::Value = serde_json::from_str(arguments_json)
            .map_err(|error| JsValue::from_str(&error.to_string()))?;
        let answer = self.inner.run(verb, &arguments).map_err(agent_error)?;
        serde_json::to_string(&answer).map_err(|error| JsValue::from_str(&error.to_string()))
    }

    /// Names of the open documents as JSON.
    pub fn documents_json(&self) -> String {
        serde_json::to_string(&self.inner.summaries()).unwrap_or_else(|_| "[]".to_owned())
    }

    /// One document as structure JSON.
    ///
    /// # Errors
    ///
    /// Returns a JavaScript object with `code` and `message` properties.
    pub fn structure_json(&self, name: Option<String>) -> Result<String, JsValue> {
        let structure = self.inner.structure(name.as_deref()).map_err(agent_error)?;
        serde_json::to_string(structure).map_err(|error| JsValue::from_str(&error.to_string()))
    }

    /// Name of the active document.
    pub fn active_name(&self) -> Option<String> {
        self.inner.active_name().map(str::to_owned)
    }
}

impl Default for AgentSession {
    fn default() -> Self {
        Self::new()
    }
}
