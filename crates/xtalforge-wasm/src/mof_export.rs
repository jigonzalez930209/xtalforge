//! MOF catalogue and builder exports (Tasks 6.11-6.17).

use wasm_bindgen::prelude::*;

use crate::common::{set, summary_object};

fn mof_error(error: xtalforge::mof::Error) -> JsValue {
    let value = js_sys::Object::new();
    let _ = js_sys::Reflect::set(
        &value,
        &JsValue::from_str("code"),
        &JsValue::from_str(error.code()),
    );
    let _ = js_sys::Reflect::set(
        &value,
        &JsValue::from_str("message"),
        &JsValue::from_str(&error.message()),
    );
    value.into()
}

/// Lists the topology catalogue.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid embedded data.
#[wasm_bindgen]
pub fn mof_topologies_info() -> Result<js_sys::Array, JsValue> {
    let array = js_sys::Array::new();
    for topology in xtalforge::mof::topologies().map_err(mof_error)? {
        let item = js_sys::Object::new();
        set(&item, "name", JsValue::from_str(&topology.name));
        set(
            &item,
            "nodes",
            JsValue::from_f64(topology.nodes.len() as f64),
        );
        set(
            &item,
            "edges",
            JsValue::from_f64(topology.edges.len() as f64),
        );
        array.push(&item.into());
    }
    Ok(array)
}

/// Lists the building blocks.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid embedded data.
#[wasm_bindgen]
pub fn mof_blocks_info() -> Result<js_sys::Array, JsValue> {
    let array = js_sys::Array::new();
    for block in xtalforge::mof::blocks().map_err(mof_error)? {
        let item = js_sys::Object::new();
        set(&item, "name", JsValue::from_str(&block.name));
        set(
            &item,
            "kind",
            JsValue::from_str(&format!("{:?}", block.kind)),
        );
        set(
            &item,
            "atoms",
            JsValue::from_f64(block.elements.len() as f64),
        );
        set(
            &item,
            "connections",
            JsValue::from_f64(block.connection_points.len() as f64),
        );
        array.push(&item.into());
    }
    Ok(array)
}

/// Builds a framework and returns its summary and warnings.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an unknown topology or block.
#[wasm_bindgen]
pub fn build_mof_info(topology: &str, node: &str, linker: &str) -> Result<JsValue, JsValue> {
    let request = xtalforge::mof::MofRequest::new(topology, node, linker);
    let built = xtalforge::mof::build_mof(&request).map_err(mof_error)?;
    let object = summary_object(&built.structure)?;
    let object: js_sys::Object = object.unchecked_into();
    let warnings = js_sys::Array::new();
    for warning in &built.warnings {
        warnings.push(&JsValue::from_str(warning));
    }
    set(&object, "warnings", warnings.into());
    Ok(object.into())
}
