//! Engine input sets (Task 18.1).
//!
//! A preset is a versioned JSON template under `data/input-sets/`. The four
//! writers only format the structure and the merged settings, so the same
//! structure, preset and overrides produce byte-identical files on every
//! surface. Overrides are merged deterministically as nested JSON objects;
//! unknown keys are rejected unless the preset marks their branch open or the
//! caller passes `"allow_unknown": true`.

use serde::{Deserialize, Serialize};
use serde_json::Value;
use xtalforge_core::Structure;

use crate::Error;

mod catalog;
mod cp2k;
mod gulp;
mod lammps;
mod merge;
mod orca;
mod qe;
mod raspa;
mod settings;
mod util;
mod vasp;

pub use catalog::{PRESET_NAMES, PresetProvenance, preset_description, preset_provenance};

/// One named input file.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct InputFile {
    /// File name as it should be written in the job directory.
    pub name: String,
    /// Complete text content of the file.
    pub text: String,
}

/// A preset identifier, for example `vasp/relax`.
pub type EnginePreset = String;

/// Overrides merged into a preset's settings.
///
/// The reserved top-level key `allow_unknown` (a boolean) permits keys that
/// the preset does not declare. Without it, an undeclared key returns
/// [`Error::UnknownOverride`].
pub type PresetOverrides = serde_json::Map<String, Value>;

/// A request for one input set.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct InputSetRequest {
    /// Preset name, as returned by [`list_presets`].
    pub preset: String,
    /// Nested settings merged over the preset defaults.
    pub overrides: PresetOverrides,
}

/// A generated input set.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct InputSetSummary {
    /// Preset that produced the files.
    pub preset: String,
    /// Files in the order declared by the preset.
    pub files: Vec<InputFile>,
}

/// A typed input-set writer.
pub trait InputSet {
    /// Builds the files for `structure`, applying `overrides`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::UnknownOverride`] for an unknown override key,
    /// [`Error::Job`] for an invalid setting value, or the error of the
    /// structure expansion.
    fn build(
        &self,
        structure: &Structure,
        overrides: &PresetOverrides,
    ) -> Result<InputSetSummary, Error>;
}

/// A concrete preset backed by one immutable catalogue entry.
#[derive(Debug, Clone)]
struct CatalogPreset {
    spec: &'static catalog::PresetSpec,
}

impl InputSet for CatalogPreset {
    fn build(
        &self,
        structure: &Structure,
        overrides: &PresetOverrides,
    ) -> Result<InputSetSummary, Error> {
        let merged = merge::merge(&self.spec.settings, overrides, &self.spec.open)?;
        let settings = settings::Settings::new(self.spec.name.clone(), merged);
        let produced = match self.spec.engine {
            catalog::Engine::Vasp => vasp::build(structure, &settings)?,
            catalog::Engine::QuantumEspresso => qe::build(structure, &settings)?,
            catalog::Engine::Cp2k => cp2k::build(structure, &settings)?,
            catalog::Engine::Lammps => lammps::build(structure, &settings)?,
            catalog::Engine::Gulp => gulp::build(structure, &settings)?,
            catalog::Engine::Raspa => raspa::build(structure, &settings)?,
            catalog::Engine::Orca => orca::build(structure, &settings)?,
        };
        Ok(InputSetSummary {
            preset: self.spec.name.clone(),
            files: assemble(&self.spec.files, produced)?,
        })
    }
}

/// Lists every available preset name in sorted order.
#[must_use]
pub fn list_presets() -> Vec<&'static str> {
    PRESET_NAMES.to_vec()
}

/// Looks up a preset by name.
///
/// The returned value is a real writer; building it never returns placeholder
/// content. `None` means the name is not registered or the embedded catalogue
/// is unreadable (which [`write_input_set`] reports as
/// [`Error::Unavailable`]).
#[must_use]
pub fn preset(name: &str) -> Option<Box<dyn InputSet>> {
    let catalogue = catalog::catalog().ok()?;
    catalogue
        .presets
        .get(name)
        .map(|spec| Box::new(CatalogPreset { spec }) as Box<dyn InputSet>)
}

/// Writes one input set from a request.
///
/// # Errors
///
/// Returns [`Error::UnknownPreset`] for an unregistered name and every error
/// of [`InputSet::build`].
pub fn write_input_set(
    structure: &Structure,
    request: &InputSetRequest,
) -> Result<InputSetSummary, Error> {
    let writer = checked_preset(&request.preset)?;
    writer.build(structure, &request.overrides)
}

fn checked_preset(name: &str) -> Result<Box<dyn InputSet>, Error> {
    catalog::catalog()
        .map_err(|error| Error::Unavailable(format!("input-set catalogue: {error}")))?;
    preset(name).ok_or_else(|| Error::UnknownPreset(name.to_owned()))
}

/// Orders `produced` as the preset declares and rejects missing or extra
/// files.
fn assemble(expected: &[String], produced: Vec<InputFile>) -> Result<Vec<InputFile>, Error> {
    let mut files = Vec::with_capacity(expected.len());
    for name in expected {
        let index = produced
            .iter()
            .position(|file| file.name == *name)
            .ok_or_else(|| Error::Missing(name.clone()))?;
        let file = produced
            .get(index)
            .cloned()
            .ok_or_else(|| Error::Missing(name.clone()))?;
        files.push(file);
    }
    if produced.len() != expected.len() {
        let extra: Vec<&str> = produced
            .iter()
            .map(|file| file.name.as_str())
            .filter(|name| !expected.iter().any(|expected| expected == name))
            .collect();
        return Err(Error::Job(format!(
            "input set produced unexpected files: {}",
            extra.join(", ")
        )));
    }
    Ok(files)
}
