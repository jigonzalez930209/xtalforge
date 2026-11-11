//! Project, parameter, document and workspace exports (Tasks 3.7-3.16).

use wasm_bindgen::prelude::*;

use crate::common::{parse_any, set, summary_object, to_js_error};

/// Returns the JSON schema of a built-in parameter preset.
///
/// # Errors
///
/// Returns a JavaScript `Error` for an unknown preset.
#[wasm_bindgen]
pub fn param_schema(name: &str) -> Result<JsValue, JsValue> {
    let set = xtalforge::core::preset(name).map_err(to_js_error)?;
    let schema = set.schema();
    let text = serde_json::to_string(&schema)
        .map_err(|error| JsValue::from_str(&format!("json_error: {error}")))?;
    js_sys::JSON::parse(&text).map_err(|error| JsValue::from_str(&format!("json_error: {error:?}")))
}

/// Packs structure text into `.xtalproj` bytes.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid text or a writer failure.
#[wasm_bindgen]
pub fn project_pack(text: &str) -> Result<Vec<u8>, JsValue> {
    let structure = parse_any(text).map_err(to_js_error)?;
    let graph = structure.bond_graph().map_err(to_js_error)?;
    let bonds = serde_json::json!(
        graph
            .bonds
            .iter()
            .map(|bond| serde_json::json!({
                "i": bond.i,
                "j": bond.j,
                "image": bond.image,
                "distance": bond.distance,
                "order": bond.order,
                "explicit": bond.explicit,
            }))
            .collect::<Vec<_>>()
    );
    let parts = xtalforge::io::project::ProjectParts::from_structure(
        &structure,
        None,
        Some(bonds),
        None,
        None,
        None,
    )
    .map_err(to_js_error)?;
    xtalforge::io::project::write_project(&parts).map_err(to_js_error)
}

/// Reads `.xtalproj` bytes and returns the structure summary.
///
/// # Errors
///
/// Returns a JavaScript `Error` for a corrupt or unsupported archive.
#[wasm_bindgen]
pub fn project_info(bytes: &[u8]) -> Result<JsValue, JsValue> {
    let parts = xtalforge::io::project::read_project(bytes).map_err(to_js_error)?;
    let structure = parts.structure().map_err(to_js_error)?;
    let object = summary_object(&structure)?;
    let object: js_sys::Object = object.unchecked_into();
    set(
        &object,
        "format_version",
        JsValue::from_f64(f64::from(xtalforge::io::project::FORMAT_VERSION)),
    );
    set(
        &object,
        "has_bonds",
        JsValue::from_bool(parts.bonds.is_some()),
    );
    set(
        &object,
        "has_sites",
        JsValue::from_bool(parts.sites.is_some()),
    );
    set(
        &object,
        "has_view",
        JsValue::from_bool(parts.view.is_some()),
    );
    set(
        &object,
        "has_session",
        JsValue::from_bool(parts.session.is_some()),
    );
    Ok(object.into())
}

/// An undoable structure document in the browser.
#[wasm_bindgen]
pub struct WasmDocument {
    inner: xtalforge::commands::Document,
}

#[wasm_bindgen]
impl WasmDocument {
    /// Creates a document from structure text.
    ///
    /// # Errors
    ///
    /// Returns a JavaScript `Error` for invalid text.
    #[wasm_bindgen(constructor)]
    pub fn new(text: &str) -> Result<WasmDocument, JsValue> {
        let structure = parse_any(text).map_err(to_js_error)?;
        Ok(Self {
            inner: xtalforge::commands::Document::new(structure),
        })
    }

    /// Moves sites by a fractional delta.
    ///
    /// # Errors
    ///
    /// Returns a JavaScript `Error` for an out-of-range site.
    pub fn move_sites(
        &mut self,
        sites: Vec<u32>,
        dx: f64,
        dy: f64,
        dz: f64,
    ) -> Result<(), JsValue> {
        let targets: Vec<(usize, [f64; 3])> = sites
            .into_iter()
            .map(|site| {
                let index = usize::try_from(site)
                    .map_err(|_| crate::common::invalid_input("bad site index"))?;
                let atom = self
                    .inner
                    .structure()
                    .sites()
                    .get(index)
                    .ok_or_else(|| crate::common::invalid_input("site out of range"))?;
                let frac = atom.frac();
                Ok((index, [frac[0] + dx, frac[1] + dy, frac[2] + dz]))
            })
            .collect::<Result<_, JsValue>>()?;
        self.inner
            .apply(
                Box::new(xtalforge::commands::MoveSites::new(targets)),
                xtalforge::commands::Origin::Script,
            )
            .map_err(to_js_error)
    }

    /// Changes one site element.
    ///
    /// # Errors
    ///
    /// Returns a JavaScript `Error` for an unknown element or site.
    pub fn set_element(&mut self, site: u32, element: &str) -> Result<(), JsValue> {
        let index = usize::try_from(site).map_err(|_| crate::common::invalid_input("bad site"))?;
        let command =
            xtalforge::commands::SetElement::new(vec![index], element).map_err(to_js_error)?;
        self.inner
            .apply(Box::new(command), xtalforge::commands::Origin::Script)
            .map_err(to_js_error)
    }

    /// Undoes the last edit and returns whether it happened.
    pub fn undo(&mut self) -> bool {
        self.inner.undo()
    }

    /// Redoes the last undone edit.
    pub fn redo(&mut self) -> bool {
        self.inner.redo()
    }

    /// Whether an undo is available.
    #[wasm_bindgen(getter)]
    pub fn can_undo(&self) -> bool {
        self.inner.can_undo()
    }

    /// Whether a redo is available.
    #[wasm_bindgen(getter)]
    pub fn can_redo(&self) -> bool {
        self.inner.can_redo()
    }

    /// Whether there are unsaved edits.
    #[wasm_bindgen(getter)]
    pub fn dirty(&self) -> bool {
        self.inner.dirty()
    }

    /// Command descriptions in the log.
    #[wasm_bindgen(getter)]
    pub fn history(&self) -> js_sys::Array {
        let array = js_sys::Array::new();
        for record in self.inner.history() {
            array.push(&JsValue::from_str(&record.description));
        }
        array
    }

    /// Writes the current structure as CIF.
    ///
    /// # Errors
    ///
    /// Returns a JavaScript `Error` for a writer failure.
    pub fn to_cif(&self) -> Result<String, JsValue> {
        xtalforge::io::write_structure(self.inner.structure(), "cif").map_err(to_js_error)
    }

    /// Returns the structure summary.
    ///
    /// # Errors
    ///
    /// Returns a JavaScript `Error` when the summary cannot be built.
    pub fn summary(&self) -> Result<JsValue, JsValue> {
        summary_object(self.inner.structure())
    }
}

/// Archives run records (a JSON array) into an in-memory ZIP.
///
/// # Errors
///
/// Returns a JavaScript `Error` for invalid records.
#[wasm_bindgen]
pub fn workspace_archive(records_json: &str) -> Result<Vec<u8>, JsValue> {
    let records: Vec<xtalforge::modules::RunRecord> = serde_json::from_str(records_json)
        .map_err(|error| crate::common::invalid_input(&format!("{error}")))?;
    xtalforge::modules::archive_records(&records).map_err(to_js_error)
}

/// Lists the run ids inside an in-memory workspace ZIP.
///
/// # Errors
///
/// Returns a JavaScript `Error` for a corrupt archive.
#[wasm_bindgen]
pub fn workspace_read_archive(bytes: &[u8]) -> Result<js_sys::Array, JsValue> {
    let records = xtalforge::modules::read_archived_records(bytes).map_err(to_js_error)?;
    let array = js_sys::Array::new();
    for record in records {
        array.push(&JsValue::from_str(&record.id));
    }
    Ok(array)
}
