//! Periodic-table reference data and element-symbol parsing.
//!
//! The data is generated from Crystal-Builder `core/elements.py`,
//! `core/pyykko.py` and `ff/vibrations.py` and embedded as a compressed JSON
//! bundle (Task 1.1). Every lookup is cached after the first decompression.

use serde::{Deserialize, Serialize};

use crate::Error;
use crate::data::{EmbeddedBundle, JsonBundle};

const ELEMENTS_SHA256: &str = "7693a55f94d0f792a9b2e17338590e33e36d41c983f0b41740ce82855d92adf6";

static ELEMENTS: JsonBundle<ElementsData> = JsonBundle::new(EmbeddedBundle::new(
    "elements",
    include_bytes!("../data/elements.json.gz"),
    ELEMENTS_SHA256,
));

/// Reference data for one element or species.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Element {
    /// Atomic number; `0` for the dummy species `X`.
    pub z: u8,
    /// Element symbol, for example `Fe`, `D` or `X`.
    pub symbol: String,
    /// Full English name.
    pub name: String,
    /// Standard atomic mass in unified atomic mass units.
    pub mass: f64,
    /// Cordero 2008 covalent radius in Angstrom.
    pub covalent_radius: f64,
    /// van der Waals radius in Angstrom (Bondi where defined).
    pub vdw_radius: f64,
    /// Default Jmol/CPK colour as 0-255 RGB.
    pub color: [u8; 3],
    /// Whether the element is a metal in the reference table.
    pub is_metal: bool,
    /// Typical valence used by the UFF atom typer, when defined.
    pub valence: Option<u8>,
    /// Pauling electronegativity, when defined.
    pub electronegativity: Option<f64>,
    /// Pyykko single-bond covalent radius in Angstrom.
    pub pyykko_single: Option<f64>,
    /// Pyykko double-bond covalent radius in Angstrom.
    pub pyykko_double: Option<f64>,
    /// Pyykko aromatic radius in Angstrom (single radius when no double one).
    pub pyykko_aromatic: Option<f64>,
}

impl Element {
    /// Whether this is the dummy connection-point species `X`.
    #[must_use]
    pub fn is_dummy(&self) -> bool {
        self.symbol == "X"
    }

    /// Whether this is an isotope species such as deuterium.
    #[must_use]
    pub fn is_isotope(&self) -> bool {
        self.symbol == "D"
    }
}

/// An isotope that the application spells with its own symbol.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Isotope {
    /// Isotope symbol, for example `D`.
    pub symbol: String,
    /// Element whose energy surface the isotope shares.
    pub element: String,
    /// Isotope mass in unified atomic mass units.
    pub mass: f64,
}

#[derive(Debug, Deserialize)]
struct ElementsData {
    schema_version: u32,
    source: String,
    elements: Vec<Element>,
    isotopes: Vec<Isotope>,
}

fn data() -> Result<&'static ElementsData, Error> {
    ELEMENTS.get()
}

/// Schema version of the embedded element bundle.
///
/// # Errors
///
/// Returns a data-bundle error if the embedded table is corrupt.
pub fn schema_version() -> Result<u32, Error> {
    Ok(data()?.schema_version)
}

/// Provenance string of the embedded element bundle.
///
/// # Errors
///
/// Returns a data-bundle error if the embedded table is corrupt.
pub fn data_source() -> Result<&'static str, Error> {
    Ok(data()?.source.as_str())
}

/// Every element and species in atomic-number order, followed by deuterium
/// and the dummy species `X`.
///
/// # Errors
///
/// Returns a data-bundle error if the embedded table is corrupt.
pub fn all() -> Result<&'static [Element], Error> {
    Ok(&data()?.elements)
}

/// Every isotope that has its own symbol.
///
/// # Errors
///
/// Returns a data-bundle error if the embedded table is corrupt.
pub fn isotopes() -> Result<&'static [Isotope], Error> {
    Ok(&data()?.isotopes)
}

/// Element by atomic number, `1..=118`; `0` returns the dummy `X`.
///
/// Deuterium is a separate species; `element_by_z(1)` returns hydrogen.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a number outside `0..=118`.
pub fn element_by_z(z: u8) -> Result<&'static Element, Error> {
    if z > 118 {
        return Err(Error::InvalidInput(format!(
            "atomic number out of range: {z}"
        )));
    }
    if z == 0 {
        return element_by_symbol("X");
    }
    all()?
        .iter()
        .find(|element| element.z == z && element.symbol != "D")
        .ok_or_else(|| Error::InvalidInput(format!("no element with atomic number {z}")))
}

/// Element by messy symbol (`Fe`, `FE`, `Fe2+`, `O_1`, `C1A`, `Ow`, `D`, `X`).
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the text has no known element symbol.
pub fn element_by_symbol(symbol: &str) -> Result<&'static Element, Error> {
    let canonical = parse_symbol(symbol)?;
    all()?
        .iter()
        .find(|element| element.symbol == canonical)
        .ok_or_else(|| Error::InvalidInput(format!("unknown element symbol: {symbol:?}")))
}

/// Canonical symbol for a `1..=118` atomic number.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a number outside `1..=118`.
pub fn symbol_from_z(z: u8) -> Result<String, Error> {
    if !(1..=118).contains(&z) {
        return Err(Error::InvalidInput(format!(
            "atomic number out of range: {z}"
        )));
    }
    Ok(element_by_z(z)?.symbol.clone())
}

/// Atomic number of a messy symbol.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the symbol is unknown.
pub fn atomic_number(symbol: &str) -> Result<u8, Error> {
    Ok(element_by_symbol(symbol)?.z)
}

/// Normalise file-shaped text into an element symbol.
///
/// The first one or two ASCII letters are tried, two letters first:
/// `Fe`, `FE`, `Fe2+`, `Fe1`, `O_1`, `C1A` and `Ow` all resolve. `Co` is
/// cobalt, while `C1` is carbon; `Kryptonite` parses as krypton because the
/// rule is prefix-based, which is why [`canonical_symbol`] exists for user
/// input.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when no known symbol is found.
pub fn parse_symbol(text: &str) -> Result<String, Error> {
    let trimmed = text.trim();
    if trimmed.is_empty() {
        return Err(Error::InvalidInput("empty element symbol".to_owned()));
    }
    let letters: String = trimmed
        .chars()
        .take_while(|character| character.is_ascii_alphabetic())
        .take(2)
        .collect();
    if letters.is_empty() {
        return Err(Error::InvalidInput(format!(
            "no element symbol in {text:?}"
        )));
    }
    if letters.len() == 2 {
        let candidate = format!(
            "{}{}",
            letters.chars().next().unwrap_or('?').to_ascii_uppercase(),
            letters.chars().nth(1).unwrap_or('?').to_ascii_lowercase()
        );
        if find_exact(&candidate)?.is_some() {
            return Ok(candidate);
        }
    }
    let one = letters
        .chars()
        .next()
        .unwrap_or('?')
        .to_ascii_uppercase()
        .to_string();
    if find_exact(&one)?.is_some() {
        return Ok(one);
    }
    Err(Error::InvalidInput(format!(
        "unknown element symbol: {text:?}"
    )))
}

/// Strict symbol reading: the whole trimmed string must be one or two letters,
/// case aside. `fe` and `FE` give `Fe`; `Fe2+`, `C1` and `Kryptonite` give
/// `None`.
///
/// # Errors
///
/// Returns a data-bundle error if the embedded table is corrupt.
pub fn canonical_symbol(text: &str) -> Result<Option<String>, Error> {
    let trimmed = text.trim();
    if trimmed.is_empty() || trimmed.len() > 2 || !trimmed.chars().all(|c| c.is_ascii_alphabetic())
    {
        return Ok(None);
    }
    let mut characters = trimmed.chars();
    let first = characters.next().unwrap_or('?').to_ascii_uppercase();
    let candidate = match characters.next() {
        Some(second) => format!("{first}{}", second.to_ascii_lowercase()),
        None => first.to_string(),
    };
    Ok(find_exact(&candidate)?.map(|_| candidate))
}

/// Whether [`canonical_symbol`] accepts the text.
///
/// # Errors
///
/// Returns a data-bundle error if the embedded table is corrupt.
pub fn is_symbol(text: &str) -> Result<bool, Error> {
    Ok(canonical_symbol(text)?.is_some())
}

/// Isotope record for a symbol, when it has one.
///
/// # Errors
///
/// Returns a data-bundle error if the embedded table is corrupt.
pub fn isotope(symbol: &str) -> Result<Option<&'static Isotope>, Error> {
    Ok(isotopes()?.iter().find(|isotope| isotope.symbol == symbol))
}

/// Mass of an isotope or element symbol in unified atomic mass units.
///
/// # Errors
///
/// Returns a data-bundle or invalid-symbol error.
pub fn mass_of(symbol: &str) -> Result<f64, Error> {
    if let Some(isotope) = isotope(parse_symbol(symbol)?.as_str())? {
        return Ok(isotope.mass);
    }
    Ok(element_by_symbol(symbol)?.mass)
}

/// Element information by atomic number (`"26"`) or messy symbol (`"Fe"`).
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the input is neither a valid number
/// nor a known symbol, or a data-bundle error if the table is corrupt.
pub fn element_info(input: &str) -> Result<Element, Error> {
    let trimmed = input.trim();
    if let Ok(z) = trimmed.parse::<u8>() {
        return Ok(element_by_z(z)?.clone());
    }
    Ok(element_by_symbol(trimmed)?.clone())
}

fn find_exact(symbol: &str) -> Result<Option<&'static Element>, Error> {
    Ok(all()?.iter().find(|element| element.symbol == symbol))
}
