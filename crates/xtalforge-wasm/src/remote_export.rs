//! Fetch-based remote providers for the browser (Task 13.3).

use wasm_bindgen::prelude::*;
use wasm_bindgen_futures::JsFuture;

async fn fetch_bytes(url: &str) -> Result<Vec<u8>, JsValue> {
    let global = js_sys::global();
    let fetch = js_sys::Reflect::get(&global, &JsValue::from_str("fetch"))?
        .dyn_into::<js_sys::Function>()
        .map_err(|_| JsValue::from_str("this runtime has no fetch"))?;
    let promise = fetch
        .call1(&JsValue::NULL, &JsValue::from_str(url))?
        .dyn_into::<js_sys::Promise>()?;
    let response = JsFuture::from(promise).await?;
    if !js_sys::Reflect::get(&response, &JsValue::from_str("ok"))?
        .as_bool()
        .unwrap_or(false)
    {
        return Err(JsValue::from_str(&format!("fetch failed: {url}")));
    }
    let array_buffer = js_sys::Reflect::get(&response, &JsValue::from_str("arrayBuffer"))?
        .dyn_into::<js_sys::Function>()?;
    let buffer = array_buffer
        .call0(&response)?
        .dyn_into::<js_sys::Promise>()?;
    let buffer = JsFuture::from(buffer).await?;
    Ok(js_sys::Uint8Array::new(&buffer).to_vec())
}

async fn fetch_text(url: &str) -> Result<String, JsValue> {
    let global = js_sys::global();
    let fetch = js_sys::Reflect::get(&global, &JsValue::from_str("fetch"))?
        .dyn_into::<js_sys::Function>()
        .map_err(|_| JsValue::from_str("this runtime has no fetch"))?;
    let promise = fetch
        .call1(&JsValue::NULL, &JsValue::from_str(url))?
        .dyn_into::<js_sys::Promise>()?;
    let response = JsFuture::from(promise).await?;
    if !js_sys::Reflect::get(&response, &JsValue::from_str("ok"))?
        .as_bool()
        .unwrap_or(false)
    {
        return Err(JsValue::from_str(&format!("fetch failed: {url}")));
    }
    let text = js_sys::Reflect::get(&response, &JsValue::from_str("text"))?
        .dyn_into::<js_sys::Function>()?;
    let body = text.call0(&response)?.dyn_into::<js_sys::Promise>()?;
    JsFuture::from(body)
        .await?
        .as_string()
        .ok_or_else(|| JsValue::from_str("remote response is not text"))
}

/// Fetches an OPTIMADE structures URL and returns its structures as JSON.
///
/// The URL must already contain the filter and pagination parameters.
#[wasm_bindgen]
pub async fn fetch_optimade(url: String) -> Result<String, JsValue> {
    let text = fetch_text(&url).await?;
    let page = xtalforge::remote::optimade_from_json(&text)
        .map_err(|error| JsValue::from_str(&error.to_string()))?;
    let documents = page
        .structures
        .iter()
        .map(|item| {
            serde_json::json!({
                "id": item.id,
                "source": item.source,
                "license": item.license,
                "structure": item.structure.to_document(),
            })
        })
        .collect::<Vec<_>>();
    serde_json::to_string(&serde_json::json!({
        "next": page.next,
        "structures": documents,
    }))
    .map_err(|error| JsValue::from_str(&error.to_string()))
}

/// Fetches a COD entry by identifier and returns its structure JSON.
#[wasm_bindgen]
pub async fn fetch_cod(id: String) -> Result<String, JsValue> {
    let url = format!("https://www.crystallography.net/cod/{id}.cif");
    let text = fetch_text(&url).await?;
    let remote = xtalforge::remote::cod_from_cif(&text, &id)
        .map_err(|error| JsValue::from_str(&error.to_string()))?;
    remote
        .structure
        .to_document()
        .to_json()
        .map_err(|error| JsValue::from_str(&error.to_string()))
}

/// Downloads a gzip data bundle on first use and returns its JSON payload.
///
/// The browser has no disk cache; the caller keeps the returned string for the
/// lifetime of the page. The SHA-256 is the digest of the **decompressed**
/// payload, the same manifest value the native loaders use.
///
/// # Errors
///
/// Returns a JavaScript `Error` when `fetch` is unavailable, the download
/// fails, the gzip payload is corrupt or the digest does not match.
#[wasm_bindgen]
pub async fn fetch_data_bundle(url: String, sha256: String) -> Result<String, JsValue> {
    let bytes = fetch_bytes(&url).await?;
    let payload = xtalforge::core::data::load_external_gzip("wasm bundle", &bytes, &sha256)
        .map_err(|error| JsValue::from_str(&error.to_string()))?;
    String::from_utf8(payload)
        .map_err(|error| JsValue::from_str(&format!("bundle is not UTF-8: {error}")))
}
