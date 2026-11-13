//! The UFF parameter table (Task 4.2.1).
//!
//! Ported from Crystal-Builder `xtal/ff/uff/params.py`, which transcribes
//! Table 1 of Rappe et al., *J. Am. Chem. Soc.* **1992**, 114, 10024, plus
//! the UFF4MOF rows of Addicoat et al., *JCTC* **2014**, 10, 880 and Coupry
//! et al., *JCTC* **2016**, 12, 5215. The rows live in `data/uff.json.gz` so a
//! published crate carries them without a large Rust literal and decompresses
//! them only on first use; the column order and units are the reference ones.
//!
//! Units: Angstrom, degrees, kcal/mol and electron charge. The public
//! calculator converts energy to eV only at the boundary.

use std::collections::BTreeMap;
use std::str::FromStr;
use std::sync::OnceLock;

use serde::{Deserialize, Serialize};

use xtalforge_core::data::{EmbeddedBundle, JsonBundle};

use crate::Error;

/// Bond-order correction coefficient, eq 3 of the paper.
pub const LAMBDA: f64 = 0.1332;
/// Bond and angle force constant prefactor, eq 6 and 13, kcal/mol/Angstrom^n.
pub const FORCE_CONSTANT: f64 = 664.12;
/// Coulomb prefactor `e^2 / (4 pi eps0)` in kcal Angstrom / mol.
pub const COULOMB: f64 = 332.0637;
/// Conversion from kcal/mol to eV.
pub const KCAL_PER_MOL_TO_EV: f64 = 0.043_364_1;
/// Conversion from kcal/mol/Angstrom^3 to GPa.
pub const KCAL_PER_MOL_A3_TO_GPA: f64 = 6.947_695_457;

const TABLE_PATH: &str = "crates/xtalforge-ff/data/uff.json.gz";
const SCHEMA_VERSION: u32 = 1;

/// Which parameters a UFF calculation may use.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub enum ParameterSet {
    /// Rappe's 1992 table only.
    #[serde(rename = "uff")]
    Uff,
    /// UFF plus the UFF4MOF and UFF4MOF-II rows.
    #[default]
    #[serde(rename = "uff4mof")]
    Uff4Mof,
}

impl ParameterSet {
    /// Stable name used by the CLI and bindings.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Uff => "uff",
            Self::Uff4Mof => "uff4mof",
        }
    }

    /// Whether a row added by UFF4MOF is allowed.
    #[must_use]
    pub fn allows(self, params: &UffParams) -> bool {
        self == Self::Uff4Mof || !params.mof_extension
    }
}

impl FromStr for ParameterSet {
    type Err = Error;

    fn from_str(text: &str) -> Result<Self, Self::Err> {
        match text.trim().to_ascii_lowercase().as_str() {
            "uff" => Ok(Self::Uff),
            "uff4mof" | "uff4mof-ii" => Ok(Self::Uff4Mof),
            other => Err(Error::InvalidInput(format!(
                "unknown UFF parameter set {other:?}; expected \"uff\" or \"uff4mof\""
            ))),
        }
    }
}

/// One row of UFF Table 1 or its UFF4MOF extension.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct UffParams {
    /// Five-character type name, for example `C_R` or `Fe6+2`.
    pub name: String,
    /// Element symbol in the IUPAC spelling (`Lr`, not the paper's `Lw`).
    pub element: String,
    /// Geometry character: `1`, `2`, `3`, `4`, `5`, `6`, `8`, `R` or empty.
    pub geometry: String,
    /// Formal oxidation state encoded in the type name, when present.
    pub oxidation_state: Option<i32>,
    /// Whether UFF4MOF fitted this row rather than Rappe's table carrying it.
    pub fitted: bool,
    /// Whether this row was added by UFF4MOF/UFF4MOF-II.
    pub mof_extension: bool,
    /// Valence bond radius in Angstrom.
    pub r1: f64,
    /// Equilibrium valence angle in degrees.
    pub theta0: f64,
    /// van der Waals distance in Angstrom.
    pub x1: f64,
    /// van der Waals well depth in kcal/mol.
    pub d1: f64,
    /// van der Waals scale, unused by the 12-6 Lennard-Jones form.
    pub zeta: f64,
    /// Effective charge `Z*`.
    pub z1: f64,
    /// sp3 torsional barrier in kcal/mol.
    pub vi: f64,
    /// sp2 torsional constant in kcal/mol.
    pub uj: f64,
    /// GMP electronegativity in eV (QEq; stored for Task 4.3).
    pub chi: f64,
    /// Idempotential / hardness in eV (QEq; stored for Task 4.3).
    pub hard: f64,
    /// QEq effective radius in Angstrom.
    pub radius: f64,
    /// `x1 / 2`, the conventional van der Waals radius in Angstrom.
    pub vdw_radius: f64,
    /// Alias of [`UffParams::d1`] in kcal/mol.
    pub vdw_energy: f64,
    /// Alias of [`UffParams::z1`].
    pub effective_charge: f64,
}

impl UffParams {
    /// Geometry character, when the name carries one.
    #[must_use]
    pub fn geometry_char(&self) -> Option<char> {
        self.geometry.chars().next()
    }

    /// Coordination number the geometry character describes, when it names
    /// one. `R` (resonant) and an empty character describe no count.
    #[must_use]
    pub fn coordination(&self) -> Option<usize> {
        match self.geometry_char()? {
            '1' => Some(2),
            '2' => Some(3),
            '3' | '4' => Some(4),
            '5' => Some(5),
            '6' => Some(6),
            '8' => Some(8),
            _ => None,
        }
    }

    /// Whether the geometry character marks a resonant (aromatic) centre.
    #[must_use]
    pub fn is_resonant(&self) -> bool {
        self.geometry_char() == Some('R')
    }

    /// Whether the type belongs to the UFF4MOF extension.
    #[must_use]
    pub fn is_mof_extension(&self) -> bool {
        self.mof_extension
    }

    /// Whether the element is a metal in the core element table.
    #[must_use]
    pub fn is_metal(&self) -> bool {
        xtalforge_core::elements::element_by_symbol(&self.element)
            .is_ok_and(|element| element.is_metal)
    }
}

#[derive(Deserialize)]
struct UffDocument {
    schema_version: u32,
    source: String,
    #[serde(default)]
    units: BTreeMap<String, String>,
    types: Vec<UffParams>,
}

#[derive(Debug)]
struct UffTable {
    schema_version: u32,
    source: String,
    units: BTreeMap<String, String>,
    params: Vec<UffParams>,
    by_name: BTreeMap<String, usize>,
    by_element: BTreeMap<String, Vec<usize>>,
}

static TABLE: OnceLock<Result<UffTable, String>> = OnceLock::new();

/// The gzip table is decompressed, checksum-verified and parsed on first use
/// (Task 15.3.2); later calls reuse the decoded document.
static DOCUMENT: JsonBundle<UffDocument> = JsonBundle::new(EmbeddedBundle::new(
    "xtalforge-ff/uff.json.gz",
    include_bytes!("../../data/uff.json.gz"),
    "593bcf48417afd9dc6f5f64ed20b469dddfa5d145ff59761c28c7cf3f18a1398",
));

fn load() -> Result<UffTable, Error> {
    let document = DOCUMENT
        .get()
        .map_err(|error| Error::Data(format!("could not read {TABLE_PATH}: {error}")))?;
    if document.schema_version != SCHEMA_VERSION {
        return Err(Error::Data(format!(
            "unsupported UFF table schema {}, expected {SCHEMA_VERSION}",
            document.schema_version
        )));
    }
    if document.types.len() != 218 {
        return Err(Error::Data(format!(
            "UFF table has {} types, expected 218",
            document.types.len()
        )));
    }
    let mut by_name = BTreeMap::new();
    let mut by_element: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (index, params) in document.types.iter().enumerate() {
        if by_name.insert(params.name.clone(), index).is_some() {
            return Err(Error::Data(format!(
                "duplicate UFF type name {}",
                params.name
            )));
        }
        by_element
            .entry(params.element.clone())
            .or_default()
            .push(index);
    }
    // The 1992 paper spells lawrencium `Lw`; the core and the IUPAC table
    // spell it `Lr`. The type name keeps the paper's spelling, but either
    // element spelling should find the rows.
    if let Some(indices) = by_element.get("Lr").cloned() {
        by_element.insert("Lw".to_owned(), indices);
    }
    Ok(UffTable {
        schema_version: document.schema_version,
        source: document.source.clone(),
        units: document.units.clone(),
        params: document.types.clone(),
        by_name,
        by_element,
    })
}

fn table() -> Result<&'static UffTable, Error> {
    match TABLE.get_or_init(|| load().map_err(|error| error.to_string())) {
        Ok(table) => Ok(table),
        Err(message) => Err(Error::Data(message.clone())),
    }
}

/// Schema version of the embedded parameter file.
///
/// # Errors
///
/// Returns [`Error::Data`] when the embedded table is corrupt.
pub fn schema_version() -> Result<u32, Error> {
    Ok(table()?.schema_version)
}

/// Provenance string of the embedded parameter file.
///
/// # Errors
///
/// Returns [`Error::Data`] when the embedded table is corrupt.
pub fn data_source() -> Result<&'static str, Error> {
    Ok(table()?.source.as_str())
}

/// Units declared by the embedded parameter file.
///
/// # Errors
///
/// Returns [`Error::Data`] when the embedded table is corrupt.
pub fn data_units() -> Result<&'static BTreeMap<String, String>, Error> {
    Ok(&table()?.units)
}

/// Every type, in table order (Rappe's 127 rows, then UFF4MOF's 91).
///
/// # Errors
///
/// Returns [`Error::Data`] when the embedded table is corrupt.
pub fn all() -> Result<&'static [UffParams], Error> {
    Ok(&table()?.params)
}

/// Parameters for a type name.
///
/// # Errors
///
/// Returns [`Error::UnknownType`] when the name is absent and [`Error::Data`]
/// when the embedded table is corrupt.
pub fn get(name: &str) -> Result<&'static UffParams, Error> {
    let table = table()?;
    let index = table
        .by_name
        .get(name)
        .ok_or_else(|| Error::UnknownType(name.to_owned()))?;
    table
        .params
        .get(*index)
        .ok_or_else(|| Error::Data(format!("UFF table index {index} is out of range")))
}

/// Types belonging to an element, in table order.
///
/// # Errors
///
/// Returns [`Error::Data`] when the embedded table is corrupt.
pub fn types_for(element: &str) -> Result<Vec<&'static UffParams>, Error> {
    let table = table()?;
    let indices = table.by_element.get(element).cloned().unwrap_or_default();
    indices
        .into_iter()
        .map(|index| {
            table
                .params
                .get(index)
                .ok_or_else(|| Error::Data(format!("UFF table index {index} is out of range")))
        })
        .collect()
}

/// Types belonging to an element under one parameter set.
///
/// # Errors
///
/// Returns the errors of [`types_for`].
pub fn types_for_set(
    element: &str,
    parameter_set: ParameterSet,
) -> Result<Vec<&'static UffParams>, Error> {
    Ok(types_for(element)?
        .into_iter()
        .filter(|params| parameter_set.allows(params))
        .collect())
}

/// Whether the table knows an element.
///
/// # Errors
///
/// Returns [`Error::Data`] when the embedded table is corrupt.
pub fn has_element(element: &str) -> Result<bool, Error> {
    Ok(table()?.by_element.contains_key(element))
}

/// All type names, in table order.
///
/// # Errors
///
/// Returns [`Error::Data`] when the embedded table is corrupt.
pub fn names() -> Result<Vec<&'static str>, Error> {
    Ok(all()?.iter().map(|params| params.name.as_str()).collect())
}
