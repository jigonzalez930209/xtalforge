//! Free solvent and counter-ion removal (Task 10.1.2).

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure, StructureBuilder};

use crate::platform::util::{bond_adjacency, components, element_counts, formula, p1_and_bonds};

#[derive(Debug, Clone, Deserialize)]
struct SolventData {
    solvents: Vec<SolventEntry>,
    cations: Vec<String>,
    anions: Vec<String>,
}

#[derive(Debug, Clone, Deserialize)]
struct SolventEntry {
    name: String,
    formula: String,
}

const SOLVENT_JSON: &str = include_str!("../../data/solvents.json");

/// Options for [`remove_solvent`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolventOptions {
    /// Remove free solvent fragments.
    pub remove_solvent: bool,
    /// Remove free counter-ions as well.
    pub remove_ions: bool,
    /// Largest non-solvent organic fragment still considered solvent.
    pub max_solvent_atoms: usize,
}

impl Default for SolventOptions {
    fn default() -> Self {
        Self {
            remove_solvent: true,
            remove_ions: false,
            max_solvent_atoms: 30,
        }
    }
}

/// One removed fragment.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RemovedFragment {
    /// Chemical formula.
    pub formula: String,
    /// Matching solvent name, when the formula is in the table.
    pub solvent: Option<String>,
    /// P1 atom indices that were removed.
    pub atoms: Vec<usize>,
    /// Whether the fragment was a counter-ion.
    pub ion: bool,
}

/// Result of solvent removal.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SolventReport {
    /// Cleaned P1 structure.
    pub structure: Structure,
    /// Fragments that were removed.
    pub removed: Vec<RemovedFragment>,
    /// Fragments that were kept because they are bound or not recognised.
    pub kept: usize,
}

/// Removes free solvent molecules and optionally free counter-ions.
///
/// The largest connected component is the framework. Every other component is
/// tested against the solvent formula table; a small organic fragment with no
/// metal is treated as solvent even when it is not in the table. Fragments
/// bonded to the framework are never removed.
///
/// # Errors
///
/// Returns the errors of the P1 expansion, bond perception or formula parsing.
pub fn remove_solvent(
    structure: &Structure,
    options: &SolventOptions,
) -> Result<SolventReport, Error> {
    let (p1, graph) = p1_and_bonds(structure)?;
    let adjacency = bond_adjacency(&graph, p1.sites().len());
    let groups = components(&adjacency);
    let main = groups
        .iter()
        .enumerate()
        .max_by_key(|(_, group)| group.len())
        .map(|(index, _)| index)
        .unwrap_or(0);
    let data: SolventData = serde_json::from_str(SOLVENT_JSON)
        .map_err(|error| Error::InvalidInput(format!("invalid solvent table: {error}")))?;
    let solvent_formulas: BTreeMap<String, String> = data
        .solvents
        .iter()
        .filter_map(|entry| {
            parse_formula(&entry.formula)
                .ok()
                .map(|counts| (formula(&counts), entry.name.clone()))
        })
        .collect();
    let ions: BTreeMap<String, ()> = data
        .cations
        .iter()
        .chain(&data.anions)
        .map(|symbol| (symbol.clone(), ()))
        .collect();
    let mut removed_atoms: Vec<usize> = Vec::new();
    let mut removed = Vec::new();
    let mut kept = 0;
    for (index, group) in groups.iter().enumerate() {
        if index == main {
            continue;
        }
        let counts = element_counts(&p1, group)?;
        let text = formula(&counts);
        let has_metal = counts.keys().any(|symbol| {
            xtalforge_core::element_info(symbol)
                .map(|info| info.is_metal)
                .unwrap_or(false)
        });
        let is_ion = group.len() == 1
            || (counts.len() == 1 && ions.contains_key(&text))
            || ions.contains_key(&text);
        let solvent_name = solvent_formulas.get(&text).cloned();
        let organic = counts.contains_key("C");
        let remove = if has_metal && !is_ion {
            false
        } else if is_ion {
            options.remove_ions
        } else {
            options.remove_solvent
                && (solvent_name.is_some() || (organic && group.len() <= options.max_solvent_atoms))
        };
        if remove {
            removed_atoms.extend(group);
            removed.push(RemovedFragment {
                formula: text,
                solvent: solvent_name,
                atoms: group.clone(),
                ion: is_ion,
            });
        } else {
            kept += 1;
        }
    }
    removed_atoms.sort_unstable();
    let mut builder = StructureBuilder::new(structure.name(), *structure.lattice());
    for (index, site) in p1.sites().iter().enumerate() {
        if removed_atoms.binary_search(&index).is_err() {
            builder = builder.site(site.clone());
        }
    }
    Ok(SolventReport {
        structure: builder
            .build()?
            .with_name(format!("{} cleaned", structure.name())),
        removed,
        kept,
    })
}

/// Parses a simple chemical formula such as `C3H7NO`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty or malformed formula.
pub fn parse_formula(text: &str) -> Result<BTreeMap<String, usize>, Error> {
    if text.trim().is_empty() {
        return Err(Error::InvalidInput(
            "a formula must not be empty".to_owned(),
        ));
    }
    let mut counts = BTreeMap::new();
    let mut symbol = String::new();
    let mut number = String::new();
    let mut flush = |symbol: &mut String, number: &mut String| -> Result<(), Error> {
        if symbol.is_empty() {
            return Ok(());
        }
        let count = if number.is_empty() {
            1
        } else {
            number
                .parse::<usize>()
                .map_err(|_| Error::InvalidInput(format!("invalid count in formula {text:?}")))?
        };
        *counts.entry(symbol.clone()).or_insert(0) += count;
        symbol.clear();
        number.clear();
        Ok(())
    };
    for character in text.chars() {
        if character.is_ascii_alphabetic() {
            if !number.is_empty() {
                flush(&mut symbol, &mut number)?;
            }
            symbol.push(character);
        } else if character.is_ascii_digit() {
            number.push(character);
        } else if character == ' ' {
            flush(&mut symbol, &mut number)?;
        } else {
            return Err(Error::InvalidInput(format!(
                "unsupported character {character:?} in formula {text:?}"
            )));
        }
    }
    flush(&mut symbol, &mut number)?;
    if counts.is_empty() {
        return Err(Error::InvalidInput(format!(
            "formula {text:?} has no elements"
        )));
    }
    Ok(counts)
}
