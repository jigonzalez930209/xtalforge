//! CoRE MOF and QMOF readers and screening (Tasks 10.6.1-10.6.3).

use std::collections::BTreeMap;
use std::path::Path;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure};
use xtalforge_symmetry::{MatcherOptions, match_structures};

use crate::platform::checker::check_mof;
use crate::platform::descriptors::framework_descriptors;
use crate::platform::mofid::mofid;
use crate::platform::porosity::{PoreOptions, pore_descriptors};

/// Database a file belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DatabaseKind {
    /// CoRE MOF 2019 or 2024.
    CoreMof,
    /// QMOF.
    Qmof,
    /// A generic directory of CIFs.
    Generic,
}

/// One database entry.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DatabaseEntry {
    /// Identifier, usually the file stem.
    pub id: String,
    /// Structure.
    pub structure: Structure,
    /// Physical properties from the database table.
    pub properties: BTreeMap<String, f64>,
    /// Source database.
    pub source: DatabaseKind,
}

/// Reader metadata.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DatabaseInfo {
    /// Source database.
    pub source: DatabaseKind,
    /// Number of structures read.
    pub entries: usize,
    /// License note.
    pub license: String,
    /// Citation.
    pub citation: String,
}

/// Reads a directory of CoRE MOF CIFs and the optional `properties.csv`.
///
/// The published archive is distributed as a directory tree; zip archives are
/// handled by the caller (extract first). This keeps the native library free
/// of a zip dependency.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the directory cannot be read or a CIF
/// does not parse, or a core error from the reader.
pub fn read_core_mof(path: &Path) -> Result<(Vec<DatabaseEntry>, DatabaseInfo), Error> {
    read_directory(path, DatabaseKind::CoreMof)
}

/// Reads a directory of QMOF CIFs.
///
/// # Errors
///
/// Returns the errors of [`read_core_mof`].
pub fn read_qmof(path: &Path) -> Result<(Vec<DatabaseEntry>, DatabaseInfo), Error> {
    read_directory(path, DatabaseKind::Qmof)
}

fn read_directory(
    path: &Path,
    source: DatabaseKind,
) -> Result<(Vec<DatabaseEntry>, DatabaseInfo), Error> {
    let mut entries = Vec::new();
    let read = std::fs::read_dir(path)
        .map_err(|error| Error::InvalidInput(format!("cannot read {}: {error}", path.display())))?;
    let mut files: Vec<std::path::PathBuf> = Vec::new();
    for item in read {
        let item = item.map_err(|error| Error::InvalidInput(error.to_string()))?;
        let file = item.path();
        if file.extension().is_some_and(|extension| extension == "cif") {
            files.push(file);
        }
    }
    files.sort();
    let properties = read_properties(path);
    for file in files {
        let text = std::fs::read_to_string(&file)
            .map_err(|error| Error::InvalidInput(error.to_string()))?;
        let structure = xtalforge_io::read_cif(&text)?;
        let id = file
            .file_stem()
            .map(|stem| stem.to_string_lossy().to_string())
            .unwrap_or_default();
        entries.push(DatabaseEntry {
            id: id.clone(),
            structure,
            properties: properties.get(&id).cloned().unwrap_or_default(),
            source,
        });
    }
    let (license, citation) = match source {
        DatabaseKind::CoreMof => (
            "CoRE MOF 2019 public v2: CC BY 4.0; cite Chung et al., J. Chem. Eng. Data 64 (2019) 5985.",
            "Chung et al. (2019), CoRE MOF 2019",
        ),
        DatabaseKind::Qmof => (
            "QMOF: CC BY 4.0; cite Rosen et al., Sci. Data 8 (2021) 217.",
            "Rosen et al. (2021), QMOF",
        ),
        DatabaseKind::Generic => ("See the source files.", "Unknown"),
    };
    let count = entries.len();
    Ok((
        entries,
        DatabaseInfo {
            source,
            entries: count,
            license: license.to_owned(),
            citation: citation.to_owned(),
        },
    ))
}

fn read_properties(path: &Path) -> BTreeMap<String, BTreeMap<String, f64>> {
    let mut output = BTreeMap::new();
    let Ok(text) = std::fs::read_to_string(path.join("properties.csv")) else {
        return output;
    };
    let mut lines = text.lines();
    let Some(header) = lines.next() else {
        return output;
    };
    let columns: Vec<&str> = header.split(',').map(str::trim).collect();
    for line in lines {
        let values: Vec<&str> = line.split(',').map(str::trim).collect();
        let Some(id) = values.first() else {
            continue;
        };
        let mut row = BTreeMap::new();
        for (column, value) in columns.iter().zip(values.iter()).skip(1) {
            if let Ok(number) = value.parse::<f64>() {
                row.insert((*column).to_owned(), number);
            }
        }
        output.insert((*id).to_owned(), row);
    }
    output
}

/// One screening record.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ScreeningRecord {
    /// Database identifier.
    pub id: String,
    /// MOFkey.
    pub mofkey: String,
    /// MOFid string.
    pub mofid: String,
    /// Framework descriptors.
    pub descriptors: crate::platform::descriptors::FrameworkDescriptors,
    /// Pore descriptors.
    pub pores: crate::platform::porosity::PoreDescriptors,
    /// Henry coefficient (dimensionless Widom average).
    pub henry: Option<f64>,
    /// Heat of adsorption in eV.
    pub heat_of_adsorption: Option<f64>,
    /// Checker findings count by severity.
    pub problems: BTreeMap<String, usize>,
}

/// Screening options.
#[derive(Debug, Clone, PartialEq)]
pub struct ScreeningOptions {
    /// Grid spacing.
    pub spacing: f64,
    /// Guest for Widom insertion.
    pub guest: String,
    /// Temperature in K.
    pub temperature: f64,
    /// Number of Widom insertions.
    pub samples: usize,
    /// Stop after this many entries (zero means all).
    pub limit: usize,
}

impl Default for ScreeningOptions {
    fn default() -> Self {
        Self {
            spacing: 0.8,
            guest: "CO2".to_owned(),
            temperature: 298.0,
            samples: 2000,
            limit: 0,
        }
    }
}

/// Cleans, deconstructs and computes descriptors and Henry coefficients for a
/// database.
///
/// The run is deterministic and writes one record per entry; per-entry errors
/// are recorded as an empty `mofkey` and the problem counts.
///
/// # Errors
///
/// Returns an error only when the options are invalid.
pub fn screen_database(
    entries: &[DatabaseEntry],
    options: &ScreeningOptions,
) -> Result<Vec<ScreeningRecord>, Error> {
    if !options.spacing.is_finite() || options.spacing <= 0.0 {
        return Err(Error::InvalidInput(
            "screening spacing must be positive and finite".to_owned(),
        ));
    }
    let limit = if options.limit == 0 {
        entries.len()
    } else {
        options.limit.min(entries.len())
    };
    let mut records = Vec::with_capacity(limit);
    for entry in entries.iter().take(limit) {
        let check = check_mof(&entry.structure).ok();
        let problems = check.map_or_else(BTreeMap::new, |report| report.counts);
        let identification = mofid(&entry.structure).ok();
        let descriptors = framework_descriptors(&entry.structure).ok();
        let pores = pore_descriptors(
            &entry.structure,
            &PoreOptions {
                spacing: options.spacing,
                ..PoreOptions::default()
            },
        )
        .ok()
        .map(|report| report.descriptors);
        let widom = crate::platform::adsorption::widom_insertion(
            &entry.structure,
            &options.guest,
            &crate::platform::adsorption::WidomOptions {
                temperature: options.temperature,
                samples: options.samples,
                seed: 0,
            },
        )
        .ok();
        let (Some(identification), Some(descriptors), Some(pores)) =
            (identification, descriptors, pores)
        else {
            records.push(ScreeningRecord {
                id: entry.id.clone(),
                mofkey: String::new(),
                mofid: String::new(),
                descriptors: empty_descriptors(&entry.structure),
                pores: empty_pores(),
                henry: widom.as_ref().map(|result| result.henry),
                heat_of_adsorption: widom.as_ref().map(|result| result.heat_of_adsorption),
                problems,
            });
            continue;
        };
        records.push(ScreeningRecord {
            id: entry.id.clone(),
            mofkey: identification.mofkey,
            mofid: identification.mofid,
            descriptors,
            pores,
            henry: widom.as_ref().map(|result| result.henry),
            heat_of_adsorption: widom.map(|result| result.heat_of_adsorption),
            problems,
        });
    }
    Ok(records)
}

fn empty_descriptors(structure: &Structure) -> crate::platform::descriptors::FrameworkDescriptors {
    crate::platform::descriptors::FrameworkDescriptors {
        formula: structure.formula(),
        n_atoms: structure.p1().map(|cell| cell.n_atoms()).unwrap_or(0),
        cell_volume: structure.lattice().volume(),
        density: structure.density().unwrap_or(0.0),
        metals: Vec::new(),
        topology: Vec::new(),
        interpenetration: 1,
        nodes: Vec::new(),
        linkers: Vec::new(),
        groups: BTreeMap::new(),
    }
}

fn empty_pores() -> crate::platform::porosity::PoreDescriptors {
    crate::platform::porosity::PoreDescriptors {
        pld: 0.0,
        lcd: 0.0,
        asa_m2_g: 0.0,
        asa_m2_cm3: 0.0,
        accessible_volume: 0.0,
        void_fraction: 0.0,
        channel_dimension: 0,
        density: 0.0,
        probe_radius: 1.2,
        spacing: 0.0,
    }
}

/// Groups structures by MOFkey, then by the structure matcher for ties.
///
/// # Errors
///
/// Returns the errors of MOFid or the matcher.
pub fn deduplicate(structures: &[Structure]) -> Result<Vec<Vec<usize>>, Error> {
    let mut keys = Vec::with_capacity(structures.len());
    for structure in structures {
        keys.push(mofid(structure)?.mofkey);
    }
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut representatives: Vec<usize> = Vec::new();
    let options = MatcherOptions {
        allow_supercell: true,
        ..MatcherOptions::default()
    };
    for (index, structure) in structures.iter().enumerate() {
        let mut found = None;
        for (group, representative) in representatives.iter().enumerate() {
            let Some(other) = structures.get(*representative) else {
                continue;
            };
            if keys.get(index) == keys.get(*representative)
                || match_structures(other, structure, &options)
                    .map(|outcome| outcome.matched)
                    .unwrap_or(false)
            {
                found = Some(group);
                break;
            }
        }
        match found {
            Some(group) => {
                if let Some(entries) = groups.get_mut(group) {
                    entries.push(index);
                }
            }
            None => {
                groups.push(vec![index]);
                representatives.push(index);
            }
        }
    }
    Ok(groups)
}
