//! Shared rebuilding helpers for the command catalogue.

use std::collections::BTreeMap;

use serde_json::Value;
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder, SymOp};

/// Rebuilds a structure with a new site list and lattice, preserving every
/// other field (name, Hall symbol, operations, bonds and metadata).
///
/// # Errors
///
/// Returns the errors of [`StructureBuilder::build`].
pub(crate) fn rebuild(
    source: &Structure,
    sites: Vec<Site>,
    lattice: Lattice,
) -> Result<Structure, Error> {
    rebuild_with(source, sites, lattice, source.hall(), source.operations())
}

/// Rebuilds a structure with new sites, lattice, Hall symbol and operations,
/// preserving name, bonds and metadata.
///
/// # Errors
///
/// Returns the errors of [`StructureBuilder::build`].
pub(crate) fn rebuild_with(
    source: &Structure,
    sites: Vec<Site>,
    lattice: Lattice,
    hall: &str,
    operations: &[SymOp],
) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new(source.name(), lattice).hall(hall);
    for operation in operations {
        builder = builder.operation(*operation);
    }
    for site in sites {
        builder = builder.site(site);
    }
    for bond in source.explicit_bonds() {
        builder = builder.explicit_bond(*bond);
    }
    for bond in source.topology_bonds() {
        builder = builder.topology_bond(*bond);
    }
    for (key, value) in source.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

/// Rebuilds a structure with a new explicit-bond list, keeping everything
/// else unchanged.
///
/// # Errors
///
/// Returns the errors of [`StructureBuilder::build`].
pub(crate) fn rebuild_explicit(
    source: &Structure,
    explicit_bonds: Vec<xtalforge_core::ExplicitBond>,
) -> Result<Structure, Error> {
    let mut builder =
        StructureBuilder::new(source.name(), *source.lattice()).hall(source.hall().to_owned());
    for operation in source.operations() {
        builder = builder.operation(*operation);
    }
    for site in source.sites() {
        builder = builder.site(site.clone());
    }
    for bond in explicit_bonds {
        builder = builder.explicit_bond(bond);
    }
    for bond in source.topology_bonds() {
        builder = builder.topology_bond(*bond);
    }
    for (key, value) in source.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

/// Returns a copy of the metadata with one entry inserted or replaced.
pub(crate) fn metadata_with(
    source: &Structure,
    key: &str,
    value: Value,
) -> BTreeMap<String, Value> {
    let mut metadata = source.metadata().clone();
    metadata.insert(key.to_owned(), value);
    metadata
}

/// Rebuilds a structure from explicit metadata.
///
/// # Errors
///
/// Returns the errors of [`StructureBuilder::build`].
pub(crate) fn rebuild_metadata(
    source: &Structure,
    sites: Vec<Site>,
    metadata: &BTreeMap<String, Value>,
) -> Result<Structure, Error> {
    let mut builder =
        StructureBuilder::new(source.name(), *source.lattice()).hall(source.hall().to_owned());
    for operation in source.operations() {
        builder = builder.operation(*operation);
    }
    for site in sites {
        builder = builder.site(site);
    }
    for bond in source.explicit_bonds() {
        builder = builder.explicit_bond(*bond);
    }
    for bond in source.topology_bonds() {
        builder = builder.topology_bond(*bond);
    }
    for (key, value) in metadata {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

/// Copies a site, replacing its element while keeping every property.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the element symbol is unknown or a
/// copied property is invalid.
pub(crate) fn retype_site(site: &Site, element: &str) -> Result<Site, Error> {
    let mut retyped =
        Site::new(site.label(), element, site.frac())?.with_occupancy(site.occupancy())?;
    if let Some(u_iso) = site.u_iso() {
        retyped = retyped.with_u_iso(u_iso)?;
    }
    if let Some(u_aniso) = site.u_aniso() {
        retyped = retyped.with_u_aniso(u_aniso)?;
    }
    if let Some(charge) = site.charge() {
        retyped = retyped.with_charge(charge)?;
    }
    if let Some(group) = site.disorder_group() {
        retyped = retyped.with_disorder_group(group);
    }
    for (key, value) in site.props() {
        retyped = retyped.with_prop(key.clone(), value.clone());
    }
    Ok(retyped)
}

/// Validates one index against a collection length.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] naming the index, the length and the kind
/// of item.
pub(crate) fn check_index(index: usize, len: usize, what: &str) -> Result<(), Error> {
    if index >= len {
        return Err(Error::InvalidInput(format!(
            "{what} index {index} is out of range for {len} items"
        )));
    }
    Ok(())
}

/// Sorts and deduplicates a list of indexes, validating every one.
///
/// # Errors
///
/// Returns the errors of [`check_index`].
pub(crate) fn checked_indices(
    indices: &[usize],
    len: usize,
    what: &str,
) -> Result<Vec<usize>, Error> {
    for index in indices {
        check_index(*index, len, what)?;
    }
    let mut sorted = indices.to_vec();
    sorted.sort_unstable();
    sorted.dedup();
    Ok(sorted)
}
