//! Runtime capabilities (Task 8.1.1).
//!
//! An assistant that quotes an engine from memory will name one that is not
//! installed here. This module answers what **this** installation can do by
//! probing it: compiled features, external binaries on `PATH`, embedded data
//! bundles and the registered formats.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_ff::{EngineKind, EngineStatus, probe};

use crate::verbs;

/// Availability of one external engine.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct EngineAvailability {
    /// Engine name as used by the force-field API.
    pub name: String,
    /// Whether the binary answered.
    pub available: bool,
    /// Path when available, or the reason it is not.
    pub detail: String,
}

/// Availability of one embedded data bundle.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct BundleAvailability {
    /// Bundle name.
    pub name: String,
    /// Whether it is embedded in this build.
    pub available: bool,
    /// What it contains.
    pub detail: String,
}

/// What this installation can do, asked at runtime.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Capabilities {
    /// Crate version.
    pub version: String,
    /// Compiled features of the agent layer.
    pub features: BTreeMap<String, bool>,
    /// External engines and whether they answer on `PATH`.
    pub engines: Vec<EngineAvailability>,
    /// Embedded data bundles.
    pub bundles: Vec<BundleAvailability>,
    /// Verbs an assistant can call.
    pub verbs: Vec<String>,
    /// Readable and writable formats.
    pub formats: Vec<FormatAvailability>,
}

/// One registered format and the directions it supports.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct FormatAvailability {
    /// Stable format name.
    pub name: String,
    /// Lower-case extensions.
    pub extensions: Vec<String>,
    /// A reader exists.
    pub readable: bool,
    /// A writer exists.
    pub writable: bool,
}

fn engine(kind: EngineKind) -> EngineAvailability {
    let status = probe(kind);
    match status {
        EngineStatus::Available { path, version } => EngineAvailability {
            name: kind.binary().to_owned(),
            available: true,
            detail: format!("{path} ({version})"),
        },
        EngineStatus::Missing { reason } => EngineAvailability {
            name: kind.binary().to_owned(),
            available: false,
            detail: reason,
        },
    }
}

fn rcsr_bundle() -> BundleAvailability {
    match xtalforge_analysis::RcsrIndex::from_embedded() {
        Ok(index) => BundleAvailability {
            name: "rcsr".to_owned(),
            available: true,
            detail: format!("{} RCSR nets, release 2019-06-01", index.len()),
        },
        Err(error) => BundleAvailability {
            name: "rcsr".to_owned(),
            available: false,
            detail: error.to_string(),
        },
    }
}

/// Probes the installation and returns its capabilities.
#[must_use]
pub fn capability_report() -> Capabilities {
    let mut features = BTreeMap::new();
    features.insert("parallel".to_owned(), cfg!(feature = "parallel"));
    features.insert("sci-form".to_owned(), cfg!(feature = "sci-form"));
    features.insert("ml".to_owned(), cfg!(feature = "ml"));
    features.insert("external".to_owned(), cfg!(feature = "external"));
    let formats = xtalforge_io::registry::FORMATS
        .iter()
        .map(|format| FormatAvailability {
            name: format.name.to_owned(),
            extensions: format
                .extensions
                .iter()
                .map(|extension| (*extension).to_owned())
                .collect(),
            readable: format.read.is_some(),
            writable: format.write.is_some(),
        })
        .collect();
    Capabilities {
        version: xtalforge_core::version().to_owned(),
        features,
        engines: vec![
            engine(EngineKind::DftbPlus),
            engine(EngineKind::Xtb),
            engine(EngineKind::Orca),
        ],
        bundles: vec![rcsr_bundle()],
        verbs: verbs::names()
            .iter()
            .map(|name| (*name).to_owned())
            .collect(),
        formats,
    }
}
