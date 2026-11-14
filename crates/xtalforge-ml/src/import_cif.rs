//! CIF directory and array import helpers (Task 21.2.3).

use std::path::Path;

use xtalforge_analysis::dataset::{Dataset, DatasetRecord};
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

use crate::import::ImportReport;

/// Imports a directory of CIF files.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unreadable directory.
pub fn import_cif_directory(path: impl AsRef<Path>, name: &str) -> Result<ImportReport, Error> {
    let mut files: Vec<std::path::PathBuf> = std::fs::read_dir(path.as_ref())
        .map_err(|error| {
            Error::InvalidInput(format!("cannot read {}: {error}", path.as_ref().display()))
        })?
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| {
            path.extension()
                .and_then(std::ffi::OsStr::to_str)
                .is_some_and(|extension| extension.eq_ignore_ascii_case("cif"))
        })
        .collect();
    files.sort();
    let mut dataset = Dataset::new(name);
    let mut skipped = 0usize;
    for file in files {
        let Ok(text) = std::fs::read_to_string(&file) else {
            skipped += 1;
            continue;
        };
        match xtalforge_io::read_structure(&text, "cif") {
            Ok(structure) => {
                let id = file
                    .file_stem()
                    .and_then(std::ffi::OsStr::to_str)
                    .unwrap_or("cif")
                    .to_owned();
                dataset.push(DatasetRecord::new(id, structure)?)?;
            }
            Err(_) => skipped += 1,
        }
    }
    let n_records = dataset.len();
    Ok(ImportReport {
        dataset,
        n_records,
        skipped,
    })
}

/// Builds a structure from explicit elements and Cartesian positions.
///
/// # Errors
///
/// Returns the builder errors.
pub fn structure_from_arrays(
    lattice: Lattice,
    elements: &[String],
    cartesian: &[[f64; 3]],
) -> Result<Structure, Error> {
    if elements.len() != cartesian.len() {
        return Err(Error::InvalidInput(
            "elements and positions must have the same length".to_owned(),
        ));
    }
    let mut builder = StructureBuilder::new("imported", lattice);
    for (index, (element, position)) in elements.iter().zip(cartesian).enumerate() {
        builder = builder.site(Site::new(
            format!("{element}{}", index + 1),
            element,
            lattice.to_frac(*position),
        )?);
    }
    builder.build()
}
