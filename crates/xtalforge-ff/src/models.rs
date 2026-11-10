//! Machine-learning model management (Task 8.4.2).
//!
//! Weights are never bundled. The catalogue records the license and homepage
//! of the three supported model families; `fetch` downloads a caller-chosen
//! artifact into the cache, records its SHA-256 and verifies it on demand.
//!
//! The download uses `curl` because the crate has no HTTP client and adding
//! one for a command that runs once per model would be a large dependency for
//! every build. A missing `curl` is a typed error with that explanation.
//!
//! The cache is `$XDG_CACHE_HOME/xtalforge/models`, or
//! `~/.cache/xtalforge/models` when the variable is unset.

use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::Error;

/// One model family the bridge supports.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub struct ModelSpec {
    /// Stable catalogue name.
    pub name: &'static str,
    /// Human-readable label.
    pub label: &'static str,
    /// License of the code and weights as published upstream.
    pub license: &'static str,
    /// Project homepage.
    pub homepage: &'static str,
    /// What the model is for.
    pub description: &'static str,
}

/// The supported model families.
pub const CATALOGUE: &[ModelSpec] = &[
    ModelSpec {
        name: "mace-mp-0",
        label: "MACE-MP-0",
        license: "MIT",
        homepage: "https://github.com/ACEsuit/mace",
        description: "Equivariant message-passing potential trained on the Materials Project trajectory",
    },
    ModelSpec {
        name: "orb-v3",
        label: "ORB v3",
        license: "Apache-2.0",
        homepage: "https://github.com/orbital-materials/orb-models",
        description: "ORB v3 universal interatomic potential",
    },
    ModelSpec {
        name: "mattersim",
        label: "MatterSim",
        license: "MIT",
        homepage: "https://github.com/microsoft/mattersim",
        description: "MatterSim deep-learning potential for materials",
    },
];

/// Metadata stored next to one downloaded artifact.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ModelRecord {
    /// Catalogue name.
    pub name: String,
    /// URL it came from.
    pub url: String,
    /// Lower-case SHA-256 of the file.
    pub sha256: String,
    /// Size in bytes.
    pub bytes: u64,
    /// License recorded at fetch time.
    pub license: String,
    /// Unix seconds when it was fetched.
    pub fetched: u64,
}

/// One catalogue entry and its cache state.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct ModelStatus {
    /// Catalogue entry.
    pub spec: ModelSpec,
    /// Metadata when the artifact is cached.
    pub record: Option<ModelRecord>,
    /// Whether a checksum verification passed on the last check; `None` when
    /// not cached.
    pub verified: Option<bool>,
    /// Path of the artifact in the cache.
    pub path: String,
}

/// Default cache directory.
#[must_use]
pub fn cache_directory() -> PathBuf {
    if let Some(path) = std::env::var_os("XTALFORGE_MODEL_CACHE")
        && !path.is_empty()
    {
        return PathBuf::from(path);
    }
    if let Some(base) = std::env::var_os("XDG_CACHE_HOME")
        && !base.is_empty()
    {
        return PathBuf::from(base).join("xtalforge").join("models");
    }
    if let Some(home) = std::env::var_os("HOME")
        && !home.is_empty()
    {
        return PathBuf::from(home)
            .join(".cache")
            .join("xtalforge")
            .join("models");
    }
    PathBuf::from("xtalforge-models")
}

fn spec(name: &str) -> Result<&'static ModelSpec, Error> {
    CATALOGUE
        .iter()
        .find(|spec| spec.name == name)
        .ok_or_else(|| Error::InvalidInput(format!("unknown model {name}")))
}

fn artifact(cache: &Path, name: &str) -> PathBuf {
    cache.join(format!("{name}.model"))
}

fn metadata(cache: &Path, name: &str) -> PathBuf {
    cache.join(format!("{name}.json"))
}

fn sha256(bytes: &[u8]) -> String {
    let digest = Sha256::digest(bytes);
    let mut out = String::with_capacity(64);
    for byte in digest {
        out.push_str(&format!("{byte:02x}"));
    }
    out
}

/// Reads the metadata of a cached model, if present.
///
/// # Errors
///
/// Returns [`Error::Data`] when the sidecar exists but is not valid JSON.
pub fn read_record(cache: &Path, name: &str) -> Result<Option<ModelRecord>, Error> {
    let path = metadata(cache, name);
    if !path.is_file() {
        return Ok(None);
    }
    let text = std::fs::read_to_string(&path)
        .map_err(|error| Error::Data(format!("{}: {error}", path.display())))?;
    let record = serde_json::from_str(&text)
        .map_err(|error| Error::Data(format!("{}: {error}", path.display())))?;
    Ok(Some(record))
}

/// Lists the catalogue with the cache state.
///
/// # Errors
///
/// Returns [`Error::Data`] for a corrupt sidecar.
pub fn list(cache: &Path) -> Result<Vec<ModelStatus>, Error> {
    let mut statuses = Vec::with_capacity(CATALOGUE.len());
    for spec in CATALOGUE {
        let record = read_record(cache, spec.name)?;
        let path = artifact(cache, spec.name);
        let verified = record.as_ref().and_then(|record| {
            std::fs::read(&path)
                .ok()
                .map(|bytes| sha256(&bytes) == record.sha256)
        });
        statuses.push(ModelStatus {
            spec: *spec,
            record,
            verified,
            path: path.display().to_string(),
        });
    }
    Ok(statuses)
}

/// Downloads one artifact and records its checksum.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown name or an empty URL, and
/// [`Error::Data`] when `curl`, the checksum or the metadata write fails.
pub fn fetch(cache: &Path, name: &str, url: &str) -> Result<ModelRecord, Error> {
    let spec = spec(name)?;
    if url.trim().is_empty() {
        return Err(Error::InvalidInput("a model URL is required".to_owned()));
    }
    std::fs::create_dir_all(cache)
        .map_err(|error| Error::Data(format!("{}: {error}", cache.display())))?;
    let partial = cache.join(format!("{name}.part"));
    let status = std::process::Command::new("curl")
        .args([
            "-L",
            "--fail",
            "--silent",
            "--show-error",
            "--output",
            &partial.display().to_string(),
            url,
        ])
        .status()
        .map_err(|error| {
            Error::Data(format!(
                "could not launch curl (install curl or download the model manually): {error}"
            ))
        })?;
    if !status.success() {
        return Err(Error::Data(format!("curl failed with {status} for {url}")));
    }
    let bytes = std::fs::read(&partial)
        .map_err(|error| Error::Data(format!("{}: {error}", partial.display())))?;
    if bytes.is_empty() {
        return Err(Error::Data(format!("{url} produced an empty file")));
    }
    let digest = sha256(&bytes);
    let target = artifact(cache, name);
    std::fs::rename(&partial, &target)
        .map_err(|error| Error::Data(format!("{}: {error}", target.display())))?;
    let fetched = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |duration| duration.as_secs());
    let record = ModelRecord {
        name: name.to_owned(),
        url: url.to_owned(),
        sha256: digest,
        bytes: bytes.len() as u64,
        license: spec.license.to_owned(),
        fetched,
    };
    write_record(cache, &record)?;
    Ok(record)
}

fn write_record(cache: &Path, record: &ModelRecord) -> Result<(), Error> {
    let path = metadata(cache, &record.name);
    let text = serde_json::to_string_pretty(record)
        .map_err(|error| Error::Data(format!("metadata: {error}")))?;
    std::fs::write(&path, text).map_err(|error| Error::Data(format!("{}: {error}", path.display())))
}

/// Verifies one cached artifact against its recorded checksum.
///
/// # Errors
///
/// Returns [`Error::Data`] when the artifact or its sidecar is missing, or the
/// checksum does not match.
pub fn verify(cache: &Path, name: &str) -> Result<ModelRecord, Error> {
    let record = read_record(cache, name)?
        .ok_or_else(|| Error::Data(format!("model {name} is not cached")))?;
    let path = artifact(cache, name);
    let bytes = std::fs::read(&path)
        .map_err(|error| Error::Data(format!("{}: {error}", path.display())))?;
    if sha256(&bytes) != record.sha256 {
        return Err(Error::Data(format!(
            "model {name} failed its SHA-256 checksum"
        )));
    }
    Ok(record)
}

/// Removes one cached artifact and its sidecar.
///
/// # Errors
///
/// Returns [`Error::Data`] when a present file cannot be removed.
pub fn remove(cache: &Path, name: &str) -> Result<bool, Error> {
    spec(name)?;
    let mut removed = false;
    for path in [artifact(cache, name), metadata(cache, name)] {
        if path.is_file() {
            std::fs::remove_file(&path)
                .map_err(|error| Error::Data(format!("{}: {error}", path.display())))?;
            removed = true;
        }
    }
    Ok(removed)
}
