//! Site edits: replacement, sorting, wrapping, adding, removing and moving.

use std::collections::BTreeMap;

use crate::edits::SplitMix64;
use crate::{Error, Lattice, Site, Structure, StructureBuilder};

/// Species replacement rules for [`replace_species`].
#[derive(Debug, Clone, PartialEq)]
pub struct SpeciesReplacement {
    /// Map from the old element symbol to the new one.
    pub map: BTreeMap<String, String>,
    /// Fraction of matching sites to replace, in `[0, 1]`.
    pub fraction: f64,
    /// Deterministic seed for the partial replacement.
    pub seed: u64,
}

/// Replaces species according to a map, optionally on a fraction of the sites.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown element, a fraction outside
/// `[0, 1]`, or a builder error.
pub fn replace_species(
    structure: &Structure,
    replacement: &SpeciesReplacement,
) -> Result<Structure, Error> {
    if !replacement.fraction.is_finite() || !(0.0..=1.0).contains(&replacement.fraction) {
        return Err(Error::InvalidInput(format!(
            "replacement fraction must be in [0, 1], got {}",
            replacement.fraction
        )));
    }
    for (from, to) in &replacement.map {
        crate::elements::element_by_symbol(from)?;
        crate::elements::element_by_symbol(to)?;
    }
    // A partial replacement is an alloy decoration of the full cell, so it
    // expands to P1 first; a full replacement keeps the asymmetric unit.
    let source = if replacement.fraction < 1.0 {
        structure.as_p1()?
    } else {
        structure.clone()
    };
    let mut indices: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (index, site) in source.sites().iter().enumerate() {
        if replacement.map.contains_key(site.element()) {
            indices
                .entry(site.element().to_owned())
                .or_default()
                .push(index);
        }
    }
    let mut selected = BTreeMap::new();
    let mut rng = SplitMix64::new(replacement.seed);
    for (element, candidates) in &indices {
        let mut shuffled = candidates.clone();
        rng.shuffle(&mut shuffled);
        let count = ((shuffled.len() as f64) * replacement.fraction).round() as usize;
        for index in shuffled.into_iter().take(count) {
            if let Some(to) = replacement.map.get(element) {
                selected.insert(index, to.clone());
            }
        }
    }
    let mut builder = StructureBuilder::new(source.name(), *source.lattice()).hall(source.hall());
    for operation in source.operations() {
        builder = builder.operation(*operation);
    }
    for (index, site) in source.sites().iter().enumerate() {
        let target = selected.get(&index).map_or(site.element(), String::as_str);
        builder = builder.site(Site::new(site.label(), target, site.frac())?);
    }
    for (key, value) in source.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

/// Site ordering key for [`sort_sites`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiteSort {
    /// Alphabetical element symbol, then coordinates.
    Element,
    /// Electronegativity (missing values last), then element.
    Electronegativity,
    /// Atomic mass, then element.
    Mass,
    /// Label, then element.
    Label,
}

/// Returns a copy with the stably sorted sites.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown element in the structure or
/// a builder error.
pub fn sort_sites(structure: &Structure, sort: SiteSort) -> Result<Structure, Error> {
    let mut sites = structure.sites().to_vec();
    let mut keyed = Vec::with_capacity(sites.len());
    for site in sites.drain(..) {
        let element = crate::elements::element_by_symbol(site.element())?;
        let key = match sort {
            SiteSort::Element => (element.z as f64, element.symbol.clone()),
            SiteSort::Electronegativity => (
                element.electronegativity.unwrap_or(f64::INFINITY),
                element.symbol.clone(),
            ),
            SiteSort::Mass => (element.mass, element.symbol.clone()),
            SiteSort::Label => (0.0, site.label().to_owned()),
        };
        keyed.push((key.0, key.1, site));
    }
    keyed.sort_by(|first, second| {
        first
            .0
            .partial_cmp(&second.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| first.1.cmp(&second.1))
    });
    let mut builder =
        StructureBuilder::new(structure.name(), *structure.lattice()).hall(structure.hall());
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for (_, _, site) in keyed {
        builder = builder.site(site);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

/// Wraps every fractional coordinate into `[0, 1)`.
///
/// # Errors
///
/// Returns a builder error.
pub fn wrap_sites(structure: &Structure) -> Result<Structure, Error> {
    let mut builder =
        StructureBuilder::new(structure.name(), *structure.lattice()).hall(structure.hall());
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for site in structure.sites() {
        builder = builder.site(site.clone().with_frac(Lattice::wrap(site.frac()))?);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

/// Adds sites to the structure.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an unknown element, or a builder error.
pub fn add_sites(structure: &Structure, sites: &[(String, [f64; 3])]) -> Result<Structure, Error> {
    let mut builder =
        StructureBuilder::new(structure.name(), *structure.lattice()).hall(structure.hall());
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for site in structure.sites() {
        builder = builder.site(site.clone());
    }
    for (element, frac) in sites {
        let label = format!("{}{}", element, structure.sites().len() + 1);
        builder = builder.site(Site::new(label, element, *frac)?);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

/// Replaces the element of one site.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an out-of-range index or unknown
/// element, or a builder error.
pub fn replace_site(
    structure: &Structure,
    index: usize,
    element: &str,
) -> Result<Structure, Error> {
    let Some(site) = structure.sites().get(index) else {
        return Err(Error::InvalidInput(format!(
            "site index {index} is out of range"
        )));
    };
    let mut builder =
        StructureBuilder::new(structure.name(), *structure.lattice()).hall(structure.hall());
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for (position, current) in structure.sites().iter().enumerate() {
        if position == index {
            builder = builder.site(Site::new(site.label(), element, site.frac())?);
        } else {
            builder = builder.site(current.clone());
        }
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

/// Moves one site to new fractional coordinates.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an out-of-range index or non-finite
/// coordinates, or a builder error.
pub fn move_site(structure: &Structure, index: usize, frac: [f64; 3]) -> Result<Structure, Error> {
    let Some(site) = structure.sites().get(index) else {
        return Err(Error::InvalidInput(format!(
            "site index {index} is out of range"
        )));
    };
    let moved = site.clone().with_frac(frac)?;
    let mut builder =
        StructureBuilder::new(structure.name(), *structure.lattice()).hall(structure.hall());
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for (position, current) in structure.sites().iter().enumerate() {
        if position == index {
            builder = builder.site(moved.clone());
        } else {
            builder = builder.site(current.clone());
        }
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

/// Removes sites by index.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an out-of-range index, or a builder
/// error.
pub fn remove_sites(structure: &Structure, indices: &[usize]) -> Result<Structure, Error> {
    let mut keep = vec![true; structure.sites().len()];
    for &index in indices {
        let Some(slot) = keep.get_mut(index) else {
            return Err(Error::InvalidInput(format!(
                "site index {index} is out of range"
            )));
        };
        *slot = false;
    }
    let mut builder =
        StructureBuilder::new(structure.name(), *structure.lattice()).hall(structure.hall());
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for (index, site) in structure.sites().iter().enumerate() {
        if keep.get(index).copied().unwrap_or(false) {
            builder = builder.site(site.clone());
        }
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}
