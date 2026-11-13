//! Remote structure providers: OPTIMADE, COD and the Materials Project
//! (Task 13.3).
//!
//! The parsing functions are available without the `net` feature so bindings
//! and tests can consume fetched JSON. The blocking HTTP client is behind the
//! `net` feature; the WASM binding uses the same parsers with a browser
//! `fetch`.

use std::path::PathBuf;

use serde_json::Value;
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

/// A structure returned by a remote provider.
#[derive(Debug, Clone, PartialEq)]
pub struct RemoteStructure {
    /// Parsed structure.
    pub structure: Structure,
    /// Provider identifier, for example `cod:1000041` or `mp-149`.
    pub id: String,
    /// Provider name.
    pub source: String,
    /// License or attribution when the provider supplies one.
    pub license: Option<String>,
}

/// One OPTIMADE page.
#[derive(Debug, Clone, PartialEq)]
pub struct OptimadePage {
    /// Structures in provider order.
    pub structures: Vec<RemoteStructure>,
    /// Opaque next-page URL, when the provider has more entries.
    pub next: Option<String>,
}

fn input(message: impl Into<String>) -> Error {
    Error::InvalidInput(message.into())
}

fn number(value: &Value, key: &str) -> Result<f64, Error> {
    value
        .as_f64()
        .ok_or_else(|| input(format!("{key} must be a number")))
}

fn vector(value: &Value, key: &str) -> Result<[f64; 3], Error> {
    let items = value
        .as_array()
        .ok_or_else(|| input(format!("{key} must be a three-number array")))?;
    if items.len() != 3 {
        return Err(input(format!("{key} must have three numbers")));
    }
    let [x, y, z] = items.as_slice() else {
        return Err(input(format!("{key} must have three numbers")));
    };
    Ok([number(x, key)?, number(y, key)?, number(z, key)?])
}

fn matrix(value: &Value, key: &str) -> Result<[[f64; 3]; 3], Error> {
    let rows = value
        .as_array()
        .ok_or_else(|| input(format!("{key} must be a 3x3 array")))?;
    if rows.len() != 3 {
        return Err(input(format!("{key} must have three rows")));
    }
    let [a, b, c] = rows.as_slice() else {
        return Err(input(format!("{key} must have three rows")));
    };
    Ok([vector(a, key)?, vector(b, key)?, vector(c, key)?])
}

/// Builds a P1 structure from rows, elements and Cartesian positions.
fn from_cartesian(
    name: &str,
    lattice: Lattice,
    elements: &[String],
    positions: &[[f64; 3]],
) -> Result<Structure, Error> {
    if elements.len() != positions.len() {
        return Err(input("remote species and position counts do not match"));
    }
    if elements.is_empty() {
        return Err(input("remote structure has no sites"));
    }
    let mut builder = StructureBuilder::new(name, lattice);
    for (index, (element, point)) in elements.iter().zip(positions).enumerate() {
        builder = builder.site(Site::new(
            format!("{element}{}", index + 1),
            element,
            lattice.to_frac(*point),
        )?);
    }
    builder.build()
}

fn optimade_entry(entry: &Value, source: &str) -> Result<RemoteStructure, Error> {
    let id = entry
        .get("id")
        .and_then(Value::as_str)
        .map(str::to_owned)
        .unwrap_or_else(|| format!("{source}-structure"));
    let attributes = entry
        .get("attributes")
        .ok_or_else(|| input("OPTIMADE entry has no attributes"))?;
    let rows = matrix(
        attributes
            .get("lattice_vectors")
            .ok_or_else(|| input("OPTIMADE entry has no lattice_vectors"))?,
        "lattice_vectors",
    )?;
    let lattice = Lattice::new(rows)?;
    let positions = attributes
        .get("cartesian_site_positions")
        .and_then(Value::as_array)
        .ok_or_else(|| input("OPTIMADE entry has no cartesian_site_positions"))?;
    let mut points = Vec::with_capacity(positions.len());
    for point in positions {
        points.push(vector(point, "cartesian_site_positions")?);
    }
    let elements = attributes
        .get("species_at_sites")
        .and_then(Value::as_array)
        .ok_or_else(|| input("OPTIMADE entry has no species_at_sites"))?
        .iter()
        .map(|item| {
            item.as_str()
                .map(str::to_owned)
                .ok_or_else(|| input("OPTIMADE species must be strings"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let structure = from_cartesian(&id, lattice, &elements, &points)?;
    let license = attributes
        .get("license")
        .and_then(Value::as_str)
        .map(str::to_owned);
    Ok(RemoteStructure {
        structure,
        id,
        source: source.to_owned(),
        license,
    })
}

/// Parses an OPTIMADE `structures` response.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for malformed JSON or missing required
/// lattice/site fields.
pub fn optimade_from_json(text: &str) -> Result<OptimadePage, Error> {
    let document: Value =
        serde_json::from_str(text).map_err(|error| input(format!("OPTIMADE JSON: {error}")))?;
    let source = document
        .get("meta")
        .and_then(|meta| meta.get("provider"))
        .and_then(Value::as_str)
        .unwrap_or("optimade");
    let data = document
        .get("data")
        .and_then(Value::as_array)
        .ok_or_else(|| input("OPTIMADE response has no data array"))?;
    let structures = data
        .iter()
        .map(|entry| optimade_entry(entry, source))
        .collect::<Result<Vec<_>, _>>()?;
    let next = document
        .get("links")
        .and_then(|links| links.get("next"))
        .and_then(Value::as_str)
        .map(str::to_owned);
    Ok(OptimadePage { structures, next })
}

/// Parses a Materials Project summary response.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for malformed JSON or a response without a
/// `structure` block.
pub fn materials_project_from_json(text: &str) -> Result<RemoteStructure, Error> {
    let document: Value = serde_json::from_str(text)
        .map_err(|error| input(format!("Materials Project JSON: {error}")))?;
    let entry = document
        .get("data")
        .and_then(Value::as_array)
        .and_then(|data| data.first())
        .ok_or_else(|| input("Materials Project response has no data"))?;
    let id = entry
        .get("material_id")
        .and_then(Value::as_str)
        .unwrap_or("mp-structure")
        .to_owned();
    let structure = entry
        .get("structure")
        .ok_or_else(|| input("Materials Project entry has no structure"))?;
    let lattice = structure
        .get("lattice")
        .ok_or_else(|| input("Materials Project structure has no lattice"))?;
    let rows = if let Some(rows_value) = lattice.get("matrix").and_then(Value::as_array) {
        matrix(&Value::Array(rows_value.clone()), "lattice.matrix")?
    } else {
        let a = number(lattice.get("a").ok_or_else(|| input("lattice.a"))?, "a")?;
        let b = number(lattice.get("b").ok_or_else(|| input("lattice.b"))?, "b")?;
        let c = number(lattice.get("c").ok_or_else(|| input("lattice.c"))?, "c")?;
        let alpha = number(
            lattice.get("alpha").ok_or_else(|| input("lattice.alpha"))?,
            "alpha",
        )?;
        let beta = number(
            lattice.get("beta").ok_or_else(|| input("lattice.beta"))?,
            "beta",
        )?;
        let gamma = number(
            lattice.get("gamma").ok_or_else(|| input("lattice.gamma"))?,
            "gamma",
        )?;
        Lattice::from_parameters(a, b, c, alpha, beta, gamma)?.matrix()
    };
    let lattice = Lattice::new(rows)?;
    let sites = structure
        .get("sites")
        .and_then(Value::as_array)
        .ok_or_else(|| input("Materials Project structure has no sites"))?;
    let mut elements = Vec::with_capacity(sites.len());
    let mut points = Vec::with_capacity(sites.len());
    for site in sites {
        elements.push(mp_species(site.get("species"))?);
        if let Some(coords) = site.get("xyz") {
            points.push(vector(coords, "structure.sites.xyz")?);
        } else if let Some(coords) = site.get("abc") {
            points.push(lattice.to_cart(vector(coords, "structure.sites.abc")?));
        } else {
            return Err(input("Materials Project site has no xyz or abc"));
        }
    }
    let parsed = from_cartesian(&id, lattice, &elements, &points)?;
    Ok(RemoteStructure {
        structure: parsed,
        id,
        source: "materials-project".to_owned(),
        license: Some("CC BY 4.0".to_owned()),
    })
}

fn mp_species(value: Option<&Value>) -> Result<String, Error> {
    match value {
        Some(Value::String(symbol)) => Ok(symbol.clone()),
        Some(Value::Object(map)) => {
            if let Some(Value::String(element)) = map.get("element") {
                return Ok(element.clone());
            }
            map.keys()
                .next()
                .cloned()
                .ok_or_else(|| input("Materials Project species object is empty"))
        }
        Some(Value::Array(items)) => {
            for item in items {
                if let Ok(symbol) = mp_species(Some(item)) {
                    return Ok(symbol);
                }
            }
            Err(input("Materials Project species array is empty"))
        }
        _ => Err(input("Materials Project site has no species")),
    }
}

/// Parses a COD CIF downloaded by `cod_fetch`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a malformed CIF.
pub fn cod_from_cif(text: &str, id: &str) -> Result<RemoteStructure, Error> {
    let structure = xtalforge_io::read_cif(text)?;
    Ok(RemoteStructure {
        structure,
        id: format!("cod:{id}"),
        source: "cod".to_owned(),
        license: Some("CC0".to_owned()),
    })
}

/// Cache directory used by the remote fetchers.
#[must_use]
pub fn cache_dir() -> PathBuf {
    if let Ok(path) = std::env::var("XTALFORGE_CACHE") {
        return PathBuf::from(path);
    }
    let home = std::env::var("HOME").unwrap_or_else(|_| ".".to_owned());
    PathBuf::from(home).join(".cache").join("xtalforge")
}

#[cfg(feature = "net")]
fn http_get(url: &str, cache_key: &str) -> Result<String, Error> {
    let path = cache_dir().join("remote").join(cache_key);
    if let Ok(text) = std::fs::read_to_string(&path) {
        return Ok(text);
    }
    let response = ureq::get(url)
        .call()
        .map_err(|error| input(format!("remote request {url} failed: {error}")))?;
    let text = response
        .into_string()
        .map_err(|error| input(format!("remote response {url} failed: {error}")))?;
    if let Some(parent) = path.parent() {
        let _ = std::fs::create_dir_all(parent);
    }
    let _ = std::fs::write(&path, &text);
    Ok(text)
}

#[cfg(feature = "net")]
fn encode(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for byte in value.bytes() {
        if byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.' | b'~') {
            out.push(char::from(byte));
        } else {
            out.push_str(&format!("%{byte:02X}"));
        }
    }
    out
}

/// Queries an OPTIMADE provider.
///
/// `base` is the provider root without `/structures`, for example
/// `https://optimade.materialsproject.org`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] on transport, HTTP or parse failures.
#[cfg(feature = "net")]
pub fn optimade_query(
    base: &str,
    filter: &str,
    limit: usize,
    offset: usize,
) -> Result<OptimadePage, Error> {
    let base = base.trim_end_matches('/');
    let url = format!(
        "{base}/structures?page_limit={limit}&page_offset={offset}&filter={}",
        encode(filter)
    );
    let key = format!("optimade-{limit}-{offset}-{:x}.json", simple_hash(&url));
    let text = http_get(&url, &key)?;
    optimade_from_json(&text)
}

/// Fetches a COD entry by numeric identifier.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] on transport, HTTP or CIF failures.
#[cfg(feature = "net")]
pub fn cod_fetch(id: &str) -> Result<RemoteStructure, Error> {
    let url = format!("https://www.crystallography.net/cod/{id}.cif");
    let text = http_get(&url, &format!("cod-{id}.cif"))?;
    cod_from_cif(&text, id)
}

/// Fetches a Materials Project structure by identifier.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] on transport, HTTP, authentication or JSON
/// failures.
#[cfg(feature = "net")]
pub fn materials_project_fetch(id: &str, api_key: &str) -> Result<RemoteStructure, Error> {
    let url = format!("https://api.materialsproject.org/materials/summary/?material_ids={id}");
    let response = ureq::get(&url)
        .set("X-API-KEY", api_key)
        .call()
        .map_err(|error| input(format!("Materials Project request failed: {error}")))?;
    let text = response
        .into_string()
        .map_err(|error| input(format!("Materials Project response failed: {error}")))?;
    let parsed = materials_project_from_json(&text)?;
    if let Some(parent) = cache_dir().join("remote").to_str() {
        let _ = std::fs::create_dir_all(parent);
    }
    Ok(parsed)
}

#[cfg(feature = "net")]
fn simple_hash(value: &str) -> u64 {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};
    let mut hasher = DefaultHasher::new();
    value.hash(&mut hasher);
    hasher.finish()
}

/// Serializes a structure as one OPTIMADE `structures` entry.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the P1 expansion fails or the
/// structure cannot be encoded.
pub fn optimade_export(structure: &Structure) -> Result<String, Error> {
    xtalforge_io::write_optimade(structure)
}
