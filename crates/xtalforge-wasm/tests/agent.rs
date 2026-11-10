//! Agent session in WebAssembly (Task 8.1).

#![cfg(target_arch = "wasm32")]
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic)]

use wasm_bindgen_test::wasm_bindgen_test;
use xtalforge_wasm::{AgentSession, agent_capabilities, agent_verbs};

const NACL: &str = include_str!("../../../tests/corpus/simple/NaCl.cif");

#[wasm_bindgen_test]
fn catalogue_and_capabilities_are_json() {
    let verbs: serde_json::Value = serde_json::from_str(&agent_verbs()).expect("verbs");
    assert!(verbs.as_array().unwrap().len() >= 20);
    let report: serde_json::Value =
        serde_json::from_str(&agent_capabilities()).expect("capabilities");
    assert!(report["version"].is_string());
}

#[wasm_bindgen_test]
fn session_loads_text_and_answers() {
    let mut session = AgentSession::new();
    let arguments = serde_json::json!({
        "text": NACL,
        "format": "cif",
        "name": "nacl"
    });
    let answer = session.run("load", &arguments.to_string()).expect("load");
    let answer: serde_json::Value = serde_json::from_str(&answer).expect("answer");
    assert_eq!(answer["kind"], "loaded");
    assert_eq!(session.active_name().as_deref(), Some("nacl"));
    let info = session
        .run("info", &serde_json::json!({}).to_string())
        .expect("info");
    let info: serde_json::Value = serde_json::from_str(&info).expect("info");
    assert_eq!(info["kind"], "info");
    assert_eq!(info["inspection"]["summary"]["formula"], "Na4 Cl4");
    let structure = session.structure_json(None).expect("structure");
    let structure: serde_json::Value = serde_json::from_str(&structure).expect("structure json");
    assert!(structure["sites"].is_array());
}

#[wasm_bindgen_test]
fn unknown_verb_returns_a_coded_error() {
    let mut session = AgentSession::new();
    let error = session.run("teleport", "{}").unwrap_err();
    let code =
        js_sys::Reflect::get(&error, &wasm_bindgen::JsValue::from_str("code")).expect("code");
    assert_eq!(code.as_string().as_deref(), Some("unknown_verb"));
}
