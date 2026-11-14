//! Embedded, versioned preset catalogue.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::OnceLock;

use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};

use super::merge;

/// Sorted names of every registered preset.
pub const PRESET_NAMES: &[&str] = &[
    "cp2k/geo-opt",
    "cp2k/mof-geo-opt",
    "cp2k/single-point",
    "gulp/lattice-relax",
    "gulp/relax",
    "lammps/uff",
    "lammps/uff4mof",
    "orca/pbe0-def2-tzvp",
    "orca/single-point",
    "qe/bands",
    "qe/nscf",
    "qe/relax",
    "qe/scf",
    "qe/vc-relax",
    "raspa/widom",
    "vasp/band",
    "vasp/dos",
    "vasp/mof",
    "vasp/relax",
    "vasp/static",
];

const SCHEMA_VERSION: u32 = 1;

const DOCUMENTS: [(&str, &str); 7] = [
    ("vasp", include_str!("../../data/input-sets/vasp.json")),
    (
        "quantum-espresso",
        include_str!("../../data/input-sets/quantum-espresso.json"),
    ),
    ("cp2k", include_str!("../../data/input-sets/cp2k.json")),
    ("lammps", include_str!("../../data/input-sets/lammps.json")),
    ("gulp", include_str!("../../data/input-sets/gulp.json")),
    ("raspa", include_str!("../../data/input-sets/raspa.json")),
    ("orca", include_str!("../../data/input-sets/orca.json")),
];

/// Provenance and template version of one preset family.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PresetProvenance {
    /// Template version, independent of the crate version.
    pub version: u32,
    /// Roadmap task that introduced the template.
    pub task: String,
    /// Human-readable source of the defaults.
    pub source: String,
    /// Documentation or reference implementation used for the template.
    pub reference: String,
    /// License of the template itself and the status of the engine data.
    pub license: String,
}

/// Catalogue engine selector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Engine {
    Vasp,
    QuantumEspresso,
    Cp2k,
    Lammps,
    Gulp,
    Raspa,
    Orca,
}

impl Engine {
    fn parse(name: &str) -> Result<Self, String> {
        match name {
            "vasp" => Ok(Self::Vasp),
            "quantum-espresso" => Ok(Self::QuantumEspresso),
            "cp2k" => Ok(Self::Cp2k),
            "lammps" => Ok(Self::Lammps),
            "gulp" => Ok(Self::Gulp),
            "raspa" => Ok(Self::Raspa),
            "orca" => Ok(Self::Orca),
            other => Err(format!("unknown engine {other:?}")),
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::Vasp => "vasp",
            Self::QuantumEspresso => "quantum-espresso",
            Self::Cp2k => "cp2k",
            Self::Lammps => "lammps",
            Self::Gulp => "gulp",
            Self::Raspa => "raspa",
            Self::Orca => "orca",
        }
    }

    fn prefix(self) -> &'static str {
        match self {
            Self::QuantumEspresso => "qe",
            other => other.name(),
        }
    }
}

/// One immutable preset.
#[derive(Debug, Clone)]
pub(crate) struct PresetSpec {
    pub(crate) name: String,
    pub(crate) engine: Engine,
    pub(crate) description: String,
    pub(crate) files: Vec<String>,
    pub(crate) settings: Map<String, Value>,
    pub(crate) open: Vec<String>,
    pub(crate) provenance: PresetProvenance,
}

/// The complete registry.
#[derive(Debug)]
pub(crate) struct Catalog {
    pub(crate) presets: BTreeMap<String, PresetSpec>,
}

static CATALOG: OnceLock<Result<Catalog, String>> = OnceLock::new();

/// Returns the parsed catalogue, loading and validating it once.
pub(crate) fn catalog() -> Result<&'static Catalog, &'static String> {
    match CATALOG.get_or_init(load) {
        Ok(catalog) => Ok(catalog),
        Err(error) => Err(error),
    }
}

/// Returns the human-readable description of a preset.
#[must_use]
pub fn preset_description(name: &str) -> Option<&'static str> {
    let catalogue = catalog().ok()?;
    catalogue
        .presets
        .get(name)
        .map(|spec| spec.description.as_str())
}

/// Returns the provenance of a preset.
#[must_use]
pub fn preset_provenance(name: &str) -> Option<PresetProvenance> {
    catalog()
        .ok()?
        .presets
        .get(name)
        .map(|spec| spec.provenance.clone())
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Document {
    schema_version: u32,
    engine: String,
    provenance: PresetProvenance,
    #[serde(default)]
    defaults: Map<String, Value>,
    #[serde(default)]
    open: Vec<String>,
    presets: BTreeMap<String, Entry>,
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Entry {
    description: String,
    files: Vec<String>,
    #[serde(default)]
    settings: Map<String, Value>,
    #[serde(default)]
    open: Vec<String>,
}

fn load() -> Result<Catalog, String> {
    let mut presets = BTreeMap::new();
    for (expected_engine, text) in DOCUMENTS {
        let document: Document = serde_json::from_str(text)
            .map_err(|error| format!("{expected_engine}: invalid JSON: {error}"))?;
        if document.schema_version != SCHEMA_VERSION {
            return Err(format!(
                "{expected_engine}: unsupported schema_version {}, expected {SCHEMA_VERSION}",
                document.schema_version
            ));
        }
        let engine = Engine::parse(&document.engine)
            .map_err(|error| format!("{expected_engine}: {error}"))?;
        if engine.name() != expected_engine {
            return Err(format!(
                "{expected_engine}: document declares engine {:?}",
                document.engine
            ));
        }
        validate_provenance(&document.provenance)
            .map_err(|error| format!("{expected_engine}: {error}"))?;
        if document.presets.is_empty() {
            return Err(format!("{expected_engine}: no presets"));
        }
        for (short_name, entry) in document.presets {
            if short_name.is_empty() || short_name.contains('/') {
                return Err(format!(
                    "{expected_engine}: invalid preset key {short_name:?}"
                ));
            }
            validate_entry(&short_name, &entry)
                .map_err(|error| format!("{expected_engine}/{short_name}: {error}"))?;
            let name = format!("{}/{short_name}", engine.prefix());
            let settings = merge::deep_merge(&document.defaults, &entry.settings)
                .into_iter()
                .collect();
            let mut open = document.open.clone();
            open.extend(entry.open.iter().cloned());
            open.sort();
            open.dedup();
            let spec = PresetSpec {
                name: name.clone(),
                engine,
                description: entry.description,
                files: entry.files,
                settings,
                open,
                provenance: document.provenance.clone(),
            };
            if presets.insert(name.clone(), spec).is_some() {
                return Err(format!("duplicate preset {name}"));
            }
        }
    }
    let expected: BTreeSet<&str> = PRESET_NAMES.iter().copied().collect();
    let actual: BTreeSet<&str> = presets.keys().map(String::as_str).collect();
    if expected != actual {
        let missing: Vec<&str> = expected.difference(&actual).copied().collect();
        let extra: Vec<&str> = actual.difference(&expected).copied().collect();
        return Err(format!(
            "preset registry mismatch: missing [{}], extra [{}]",
            missing.join(", "),
            extra.join(", ")
        ));
    }
    Ok(Catalog { presets })
}

fn validate_provenance(provenance: &PresetProvenance) -> Result<(), String> {
    if provenance.version == 0 {
        return Err("provenance.version must be at least 1".to_owned());
    }
    for (name, value) in [
        ("task", &provenance.task),
        ("source", &provenance.source),
        ("reference", &provenance.reference),
        ("license", &provenance.license),
    ] {
        if value.trim().is_empty() {
            return Err(format!("provenance.{name} must not be empty"));
        }
    }
    Ok(())
}

fn validate_entry(short_name: &str, entry: &Entry) -> Result<(), String> {
    if entry.description.trim().is_empty() {
        return Err(format!("{short_name}: description must not be empty"));
    }
    if entry.files.is_empty() {
        return Err(format!("{short_name}: files must not be empty"));
    }
    let mut names = BTreeSet::new();
    for file in &entry.files {
        let valid = !file.is_empty()
            && file != "."
            && file != ".."
            && !file.contains('/')
            && !file.contains('\\')
            && !file.chars().any(char::is_control);
        if !valid {
            return Err(format!("{short_name}: invalid file name {file:?}"));
        }
        if !names.insert(file.as_str()) {
            return Err(format!("{short_name}: duplicate file name {file:?}"));
        }
    }
    Ok(())
}
