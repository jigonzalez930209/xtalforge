//! Fragment library (Task 6.2.1).
//!
//! Ported from Crystal-Builder `xtal/build/data/fragments.json` (MIT): a
//! molecule worth not typing twice, with a name, SMILES, category and
//! description. The geometry is built by [`crate::parse_smiles`] until the 3D
//! embedder lands.

use serde::Deserialize;
use xtalforge_core::Error;

#[derive(Debug, Deserialize)]
struct Data {
    fragments: Vec<Fragment>,
}

/// One library entry.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct Fragment {
    /// Display name.
    pub name: String,
    /// SMILES string.
    pub smiles: String,
    /// Category.
    pub category: String,
    /// Description.
    pub description: String,
}

/// The embedded fragment library.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] if the embedded data is invalid.
pub fn fragments() -> Result<Vec<Fragment>, Error> {
    let data: Data = serde_json::from_str(include_str!("../data/fragments.json"))
        .map_err(|error| Error::InvalidInput(format!("fragment data: {error}")))?;
    Ok(data.fragments)
}

/// Finds an entry by case-insensitive name.
///
/// # Errors
///
/// Returns the errors of [`fragments`].
pub fn find_fragment(name: &str) -> Result<Option<Fragment>, Error> {
    Ok(fragments()?
        .into_iter()
        .find(|fragment| fragment.name.eq_ignore_ascii_case(name)))
}
