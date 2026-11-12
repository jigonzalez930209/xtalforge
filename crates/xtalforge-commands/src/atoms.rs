//! Undoable atom edits on the asymmetric unit (Task 3.10.1).

use xtalforge_core::{
    Error, ExplicitBond, Site, Structure, StructureBuilder, TopologyBond,
    elements::element_by_symbol,
};

use crate::command::Command;
use crate::edit::{checked_indices, rebuild, retype_site};

/// Appends sites to the asymmetric unit.
#[derive(Debug, Clone)]
pub struct AddSites {
    sites: Vec<Site>,
    description: String,
}

impl AddSites {
    /// Creates the command from fully built sites.
    #[must_use]
    pub fn new(sites: Vec<Site>) -> Self {
        let description = match sites.as_slice() {
            [] => "Add sites".to_owned(),
            [site] => format!("Add {}", site.element()),
            _ => format!("Add {} sites", sites.len()),
        };
        Self { sites, description }
    }

    /// Sites the command appends.
    #[must_use]
    pub fn sites(&self) -> &[Site] {
        &self.sites
    }
}

impl Command for AddSites {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        let mut sites = structure.sites().to_vec();
        sites.extend(self.sites.iter().cloned());
        rebuild(structure, sites, *structure.lattice())
    }

    fn description(&self) -> String {
        self.description.clone()
    }
}

/// Removes sites from the asymmetric unit.
#[derive(Debug, Clone)]
pub struct DeleteSites {
    indices: Vec<usize>,
    description: String,
}

impl DeleteSites {
    /// Creates the command from asymmetric-unit indexes.
    #[must_use]
    pub fn new(indices: Vec<usize>) -> Self {
        let description = if indices.len() == 1 {
            "Delete site".to_owned()
        } else {
            format!("Delete {} sites", indices.len())
        };
        Self {
            indices,
            description,
        }
    }

    /// Indexes the command removes.
    #[must_use]
    pub fn indices(&self) -> &[usize] {
        &self.indices
    }
}

impl Command for DeleteSites {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        let removed = checked_indices(&self.indices, structure.sites().len(), "site")?;
        if removed.is_empty() {
            return Ok(structure.clone());
        }
        let old_len = structure.sites().len();
        let mut mapping = vec![None; old_len];
        let mut next = 0;
        for old in 0..old_len {
            if removed.binary_search(&old).is_ok() {
                continue;
            }
            if let Some(slot) = mapping.get_mut(old) {
                *slot = Some(next);
            }
            next += 1;
        }
        let sites: Vec<Site> = structure
            .sites()
            .iter()
            .enumerate()
            .filter(|(index, _)| removed.binary_search(index).is_err())
            .map(|(_, site)| site.clone())
            .collect();
        remapped_structure(structure, sites, &mapping)
    }

    fn description(&self) -> String {
        self.description.clone()
    }
}

/// Changes the element of one or more sites, keeping their properties.
#[derive(Debug, Clone)]
pub struct SetElement {
    indices: Vec<usize>,
    element: String,
    description: String,
}

impl SetElement {
    /// Creates the command.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the element symbol is unknown.
    pub fn new(indices: Vec<usize>, element: &str) -> Result<Self, Error> {
        let canonical = element_by_symbol(element)?.symbol.clone();
        let description = if indices.len() == 1 {
            format!("Change element to {canonical}")
        } else {
            format!("Change {} sites to {canonical}", indices.len())
        };
        Ok(Self {
            indices,
            element: canonical,
            description,
        })
    }

    /// Indexes the command edits.
    #[must_use]
    pub fn indices(&self) -> &[usize] {
        &self.indices
    }

    /// Canonical element symbol.
    #[must_use]
    pub fn element(&self) -> &str {
        &self.element
    }
}

impl Command for SetElement {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        let indices = checked_indices(&self.indices, structure.sites().len(), "site")?;
        let mut sites = structure.sites().to_vec();
        for index in indices {
            let Some(site) = structure.sites().get(index) else {
                continue;
            };
            let retyped = retype_site(site, &self.element)?;
            if let Some(slot) = sites.get_mut(index) {
                *slot = retyped;
            }
        }
        rebuild(structure, sites, *structure.lattice())
    }

    fn description(&self) -> String {
        self.description.clone()
    }
}

/// Moves sites to absolute fractional-coordinate targets.
#[derive(Debug, Clone)]
pub struct MoveSites {
    targets: Vec<(usize, [f64; 3])>,
    description: String,
}

impl MoveSites {
    /// Creates the command from `(site index, fractional target)` pairs.
    #[must_use]
    pub fn new(targets: Vec<(usize, [f64; 3])>) -> Self {
        let description = if targets.len() == 1 {
            "Move site".to_owned()
        } else {
            format!("Move {} sites", targets.len())
        };
        Self {
            targets,
            description,
        }
    }

    /// Targets the command applies.
    #[must_use]
    pub fn targets(&self) -> &[(usize, [f64; 3])] {
        &self.targets
    }
}

impl Command for MoveSites {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        let mut sites = structure.sites().to_vec();
        for (index, target) in &self.targets {
            if *index >= structure.sites().len() {
                return Err(Error::InvalidInput(format!(
                    "site index {index} is out of range for {} items",
                    structure.sites().len()
                )));
            }
            if !target.iter().all(|value| value.is_finite()) {
                return Err(Error::InvalidInput(
                    "site coordinates must be finite".to_owned(),
                ));
            }
            let Some(site) = structure.sites().get(*index) else {
                continue;
            };
            let moved = site.clone().with_frac(*target)?;
            if let Some(slot) = sites.get_mut(*index) {
                *slot = moved;
            }
        }
        rebuild(structure, sites, *structure.lattice())
    }

    fn description(&self) -> String {
        self.description.clone()
    }
}

/// Changes the occupancy of one site.
#[derive(Debug, Clone)]
pub struct SetOccupancy {
    index: usize,
    occupancy: f64,
}

impl SetOccupancy {
    /// Creates the command.
    #[must_use]
    pub fn new(index: usize, occupancy: f64) -> Self {
        Self { index, occupancy }
    }

    /// Site the command edits.
    #[must_use]
    pub fn index(&self) -> usize {
        self.index
    }

    /// Requested occupancy.
    #[must_use]
    pub fn occupancy(&self) -> f64 {
        self.occupancy
    }
}

impl Command for SetOccupancy {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        if self.index >= structure.sites().len() {
            return Err(Error::InvalidInput(format!(
                "site index {} is out of range for {} items",
                self.index,
                structure.sites().len()
            )));
        }
        let mut sites = structure.sites().to_vec();
        let Some(site) = structure.sites().get(self.index) else {
            return Ok(structure.clone());
        };
        let updated = site.clone().with_occupancy(self.occupancy)?;
        if let Some(slot) = sites.get_mut(self.index) {
            *slot = updated;
        }
        rebuild(structure, sites, *structure.lattice())
    }

    fn description(&self) -> String {
        format!("Set occupancy {}", self.occupancy)
    }
}

/// Builds the structure after a removal, remapping explicit bonds.
fn remapped_structure(
    structure: &Structure,
    sites: Vec<Site>,
    mapping: &[Option<usize>],
) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new(structure.name(), *structure.lattice())
        .hall(structure.hall().to_owned());
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for site in sites {
        builder = builder.site(site);
    }
    for bond in structure.explicit_bonds() {
        if let (Some(i), Some(j)) = (
            mapping.get(bond.i).copied().flatten(),
            mapping.get(bond.j).copied().flatten(),
        ) {
            builder = builder.explicit_bond(ExplicitBond { i, j, ..*bond });
        }
    }
    for bond in structure.topology_bonds() {
        if let (Some(i), Some(j)) = (
            mapping.get(bond.i).copied().flatten(),
            mapping.get(bond.j).copied().flatten(),
        ) {
            builder = builder.topology_bond(TopologyBond { i, j, ..*bond });
        }
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

/// Replaces species by a map, optionally on a fraction of the sites.
#[derive(Debug, Clone)]
pub struct ReplaceSpecies {
    replacement: xtalforge_core::SpeciesReplacement,
    description: String,
}

impl ReplaceSpecies {
    /// Creates the command.
    #[must_use]
    pub fn new(replacement: xtalforge_core::SpeciesReplacement) -> Self {
        let description = if replacement.fraction >= 1.0 {
            "Replace species".to_owned()
        } else {
            format!(
                "Replace species on {:.0}% of the sites",
                replacement.fraction * 100.0
            )
        };
        Self {
            replacement,
            description,
        }
    }

    /// Replacement rules.
    #[must_use]
    pub fn replacement(&self) -> &xtalforge_core::SpeciesReplacement {
        &self.replacement
    }
}

impl Command for ReplaceSpecies {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        xtalforge_core::replace_species(structure, &self.replacement)
    }

    fn description(&self) -> String {
        self.description.clone()
    }
}

/// Sorts the sites deterministically.
#[derive(Debug, Clone)]
pub struct SortSites {
    sort: xtalforge_core::SiteSort,
}

impl SortSites {
    /// Creates the command.
    #[must_use]
    pub fn new(sort: xtalforge_core::SiteSort) -> Self {
        Self { sort }
    }

    /// Ordering key.
    #[must_use]
    pub fn sort(&self) -> xtalforge_core::SiteSort {
        self.sort
    }
}

impl Command for SortSites {
    fn apply(&self, structure: &Structure) -> Result<Structure, Error> {
        xtalforge_core::sort_sites(structure, self.sort)
    }

    fn description(&self) -> String {
        "Sort sites".to_owned()
    }
}
