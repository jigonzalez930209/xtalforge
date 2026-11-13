//! Project files: a ZIP of small readable JSON parts (Task 3.7).
//!
//! The container is a ZIP archive with deflate compression. Parts:
//! `structure.cif`, `cell.json`, `bonds.json`, `sites.json`, `view.json` and
//! `session.json`, plus `format.json` carrying the independent
//! [`FORMAT_VERSION`].
//!
//! Readers accept every older format version through the explicit migration
//! below; a newer version is refused instead of being misread.

use std::io::{Cursor, Read, Write};

use serde_json::Value;
use xtalforge_core::{Error, Provenance, Structure};
use zip::write::SimpleFileOptions;

use crate::cif::{read_cif, write_cif};
use crate::provenance::PROJECT_PROVENANCE_KEY;

/// Current project format version, independent of the crate version.
pub const FORMAT_VERSION: u32 = 1;

/// The readable parts of a project.
#[derive(Debug, Clone, PartialEq)]
pub struct ProjectParts {
    /// The structure as CIF.
    pub structure_cif: String,
    /// Bare cell when the project has no atoms.
    pub cell: Option<Value>,
    /// Perceived and explicit bonds.
    pub bonds: Option<Value>,
    /// Per-site properties a CIF cannot hold.
    pub sites: Option<Value>,
    /// View settings.
    pub view: Option<Value>,
    /// Selection, measurements and provenance.
    pub session: Option<Value>,
}

impl ProjectParts {
    /// Builds a project from a structure and optional parts.
    ///
    /// # Errors
    ///
    /// Returns the errors of the CIF writer.
    pub fn from_structure(
        structure: &Structure,
        cell: Option<Value>,
        bonds: Option<Value>,
        sites: Option<Value>,
        view: Option<Value>,
        session: Option<Value>,
    ) -> Result<Self, Error> {
        Ok(Self {
            structure_cif: write_cif(structure)?,
            cell,
            bonds,
            sites,
            view,
            session,
        })
    }

    /// Parses the structure part.
    ///
    /// # Errors
    ///
    /// Returns the errors of the CIF reader.
    pub fn structure(&self) -> Result<Structure, Error> {
        read_cif(&self.structure_cif)
    }

    /// Stores a provenance record in the `provenance` key of `session.json`.
    ///
    /// An existing record is replaced and every other session entry is kept.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the record is invalid, the session
    /// part is present but not a JSON object, or the record cannot be
    /// serialized.
    pub fn with_provenance(mut self, record: &Provenance) -> Result<Self, Error> {
        record.validate()?;
        let value = serde_json::to_value(record).map_err(|error| {
            Error::InvalidInput(format!("provenance: could not serialize record: {error}"))
        })?;
        let mut session = match self.session.take() {
            Some(Value::Object(session)) => session,
            Some(_) => {
                return Err(Error::InvalidInput(
                    "project session part is not a JSON object".to_owned(),
                ));
            }
            None => serde_json::Map::new(),
        };
        session.insert(PROJECT_PROVENANCE_KEY.to_owned(), value);
        self.session = Some(Value::Object(session));
        Ok(self)
    }

    /// Reads the provenance record stored in `session.json`, if any.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the record is present but invalid.
    pub fn provenance(&self) -> Result<Option<Provenance>, Error> {
        let Some(session) = &self.session else {
            return Ok(None);
        };
        let Some(value) = session.get(PROJECT_PROVENANCE_KEY) else {
            return Ok(None);
        };
        Provenance::from_value(value).map(Some)
    }
}

/// Serializes a project to ZIP bytes.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when a part cannot be serialized or the
/// archive cannot be written.
pub fn write_project(parts: &ProjectParts) -> Result<Vec<u8>, Error> {
    let mut buffer = Cursor::new(Vec::new());
    {
        let mut archive = zip::ZipWriter::new(&mut buffer);
        let options = SimpleFileOptions::default()
            .compression_method(zip::CompressionMethod::Deflated)
            .unix_permissions(0o644);
        let entries = [
            ("structure.cif", parts.structure_cif.clone().into_bytes()),
            (
                "format.json",
                format!("{{\"format_version\":{FORMAT_VERSION}}}\n").into_bytes(),
            ),
        ];
        for (name, data) in entries {
            archive.start_file(name, options).map_err(zip_error)?;
            archive.write_all(&data).map_err(zip_error)?;
        }
        for (name, part) in [
            ("cell.json", &parts.cell),
            ("bonds.json", &parts.bonds),
            ("sites.json", &parts.sites),
            ("view.json", &parts.view),
            ("session.json", &parts.session),
        ] {
            if let Some(value) = part {
                let text = serde_json::to_vec_pretty(value).map_err(|error| {
                    Error::InvalidInput(format!("could not serialize {name}: {error}"))
                })?;
                archive.start_file(name, options).map_err(zip_error)?;
                archive.write_all(&text).map_err(zip_error)?;
            }
        }
        archive.finish().map_err(zip_error)?;
    }
    Ok(buffer.into_inner())
}

/// Reads a project from ZIP bytes.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a corrupt archive, a missing structure
/// part, an unsupported version or an invalid JSON part.
pub fn read_project(bytes: &[u8]) -> Result<ProjectParts, Error> {
    let cursor = Cursor::new(bytes);
    let mut archive = zip::ZipArchive::new(cursor).map_err(zip_error)?;
    let mut entries: std::collections::BTreeMap<String, Vec<u8>> =
        std::collections::BTreeMap::new();
    for index in 0..archive.len() {
        let mut file = archive.by_index(index).map_err(zip_error)?;
        let Some(name) = file.enclosed_name() else {
            return Err(Error::InvalidInput(
                "project part name escapes the archive".to_owned(),
            ));
        };
        let name = name.to_string_lossy().into_owned();
        let mut data = Vec::new();
        file.read_to_end(&mut data).map_err(zip_error)?;
        entries.insert(name, data);
    }
    let structure_cif = entries
        .get("structure.cif")
        .ok_or_else(|| Error::InvalidInput("project has no structure.cif".to_owned()))?;
    let structure_cif = String::from_utf8(structure_cif.clone())
        .map_err(|_| Error::InvalidInput("structure.cif is not UTF-8".to_owned()))?;
    let version = entries
        .get("format.json")
        .and_then(|text| std::str::from_utf8(text).ok())
        .and_then(|text| serde_json::from_str::<Value>(text).ok())
        .and_then(|value| {
            value
                .get("format_version")
                .and_then(serde_json::Value::as_u64)
        })
        .unwrap_or(1);
    if version > u64::from(FORMAT_VERSION) {
        return Err(Error::InvalidInput(format!(
            "project format_version {version} is newer than {FORMAT_VERSION}"
        )));
    }
    let part = |name: &str| -> Result<Option<Value>, Error> {
        entries
            .get(name)
            .map(|text| {
                serde_json::from_slice(text)
                    .map_err(|error| Error::InvalidInput(format!("invalid {name}: {error}")))
            })
            .transpose()
    };
    Ok(ProjectParts {
        structure_cif,
        cell: part("cell.json")?,
        bonds: part("bonds.json")?,
        sites: part("sites.json")?,
        view: part("view.json")?,
        session: part("session.json")?,
    })
}

fn zip_error(error: impl std::fmt::Display) -> Error {
    Error::InvalidInput(format!("project archive: {error}"))
}
