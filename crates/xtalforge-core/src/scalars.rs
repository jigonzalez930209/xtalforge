//! Named per-atom and per-bond scalars (Task 1.14.2).

use crate::elements::element_by_symbol;
use crate::{Error, Structure};

/// Bond lengths in the perceived graph.
///
/// # Errors
///
/// Returns the errors of the P1 expansion or bond perception.
pub fn bond_lengths(structure: &Structure) -> Result<Vec<f64>, Error> {
    Ok(structure
        .bond_graph()?
        .bonds
        .iter()
        .map(|bond| bond.distance)
        .collect())
}

/// Coordination number of every P1 atom.
///
/// # Errors
///
/// Returns the errors of the P1 expansion or bond perception.
pub fn coordination_numbers(structure: &Structure) -> Result<Vec<usize>, Error> {
    let graph = structure.bond_graph()?;
    Ok((0..graph.n_atoms())
        .map(|index| graph.degree(index))
        .collect())
}

/// Formal charge of every asymmetric site, when present.
#[must_use]
pub fn atom_charges(structure: &Structure) -> Vec<Option<f64>> {
    structure.sites().iter().map(|site| site.charge()).collect()
}

/// Molar mass of every asymmetric site.
///
/// # Errors
///
/// Returns an error if an element is missing from the table.
pub fn atom_masses(structure: &Structure) -> Result<Vec<f64>, Error> {
    structure
        .sites()
        .iter()
        .map(|site| Ok(element_by_symbol(site.element())?.mass))
        .collect()
}

/// Metadata for one named scalar.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ScalarInfo {
    /// Stable registry name.
    pub name: String,
    /// One-line description.
    pub description: String,
    /// Unit, or an empty string for a dimensionless value.
    pub unit: String,
}

/// A function that computes one scalar for every atom or bond.
pub type ScalarFn = fn(&Structure) -> Result<Vec<f64>, Error>;

struct ScalarEntry {
    info: ScalarInfo,
    compute: ScalarFn,
}

/// A registry of named scalars (Task 8.3.1).
///
/// The crate ships the four built-in scalars; a plugin can add its own and a
/// front end can list them. Registration is explicit: there is no global
/// mutable registry, so tests and servers do not interfere.
#[derive(Default)]
pub struct ScalarRegistry {
    entries: std::collections::BTreeMap<String, ScalarEntry>,
}

impl std::fmt::Debug for ScalarRegistry {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("ScalarRegistry")
            .field("names", &self.names())
            .finish()
    }
}

fn charges(structure: &Structure) -> Result<Vec<f64>, Error> {
    Ok(atom_charges(structure)
        .into_iter()
        .map(|charge| charge.unwrap_or(f64::NAN))
        .collect())
}

fn coordination(structure: &Structure) -> Result<Vec<f64>, Error> {
    Ok(coordination_numbers(structure)?
        .into_iter()
        .map(|number| number as f64)
        .collect())
}

fn builtin(name: &str, description: &str, unit: &str, compute: ScalarFn) -> ScalarEntry {
    ScalarEntry {
        info: ScalarInfo {
            name: name.to_owned(),
            description: description.to_owned(),
            unit: unit.to_owned(),
        },
        compute,
    }
}

impl ScalarRegistry {
    /// A registry with the four built-in scalars.
    #[must_use]
    pub fn new() -> Self {
        let mut registry = Self::default();
        for entry in [
            builtin(
                "bond_lengths",
                "Perceived bond lengths",
                "Angstrom",
                bond_lengths,
            ),
            builtin(
                "coordination_numbers",
                "Perceived coordination number per P1 atom",
                "",
                coordination,
            ),
            builtin(
                "atom_charges",
                "Stored charge per asymmetric site; NaN when absent",
                "e",
                charges,
            ),
            builtin(
                "atom_masses",
                "Standard atomic mass per asymmetric site",
                "u",
                atom_masses,
            ),
        ] {
            registry.entries.insert(entry.info.name.clone(), entry);
        }
        registry
    }

    /// An empty registry.
    #[must_use]
    pub fn empty() -> Self {
        Self::default()
    }

    /// Registers a scalar.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an empty or duplicate name.
    pub fn register(&mut self, info: ScalarInfo, compute: ScalarFn) -> Result<(), Error> {
        if info.name.trim().is_empty() {
            return Err(Error::InvalidInput(
                "a scalar needs a non-empty name".to_owned(),
            ));
        }
        if self.entries.contains_key(&info.name) {
            return Err(Error::InvalidInput(format!(
                "scalar {} is already registered",
                info.name
            )));
        }
        self.entries
            .insert(info.name.clone(), ScalarEntry { info, compute });
        Ok(())
    }

    /// Computes one named scalar.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an unknown name or the scalar's own
    /// error.
    pub fn compute(&self, name: &str, structure: &Structure) -> Result<Vec<f64>, Error> {
        let entry = self
            .entries
            .get(name)
            .ok_or_else(|| Error::InvalidInput(format!("unknown scalar {name}")))?;
        (entry.compute)(structure)
    }

    /// Metadata for one scalar.
    #[must_use]
    pub fn info(&self, name: &str) -> Option<&ScalarInfo> {
        self.entries.get(name).map(|entry| &entry.info)
    }

    /// Registered names, sorted.
    #[must_use]
    pub fn names(&self) -> Vec<&str> {
        self.entries.keys().map(String::as_str).collect()
    }
}
