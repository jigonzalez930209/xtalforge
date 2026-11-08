//! Structure and site helpers shared by the force-field modules.
//!
//! Rebuilding a structure is the one operation charges and optimizers both
//! need: replace the sites and/or the lattice while keeping the name, Hall
//! symbol, symmetry operations, bond overrides and metadata.  It lives here
//! so there is one copy of that rule.

use xtalforge_core::{Lattice, Site, Structure, StructureBuilder};

use crate::Error;

/// Rebuilds `structure` with new sites and lattice, preserving everything
/// else.
pub(crate) fn rebuild(
    structure: &Structure,
    lattice: Lattice,
    sites: Vec<Site>,
) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new(structure.name(), lattice).hall(structure.hall());
    for site in sites {
        builder = builder.site(site);
    }
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for bond in structure.explicit_bonds() {
        builder = builder.explicit_bond(*bond);
    }
    for bond in structure.topology_bonds() {
        builder = builder.topology_bond(*bond);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build().map_err(Error::from)
}

/// Copy of `structure` with new fractional coordinates for its sites and a
/// new lattice.
pub(crate) fn with_frac(
    structure: &Structure,
    lattice: Lattice,
    frac: &[[f64; 3]],
) -> Result<Structure, Error> {
    if frac.len() != structure.sites().len() {
        return Err(Error::InvalidInput(format!(
            "expected {} site coordinates, got {}",
            structure.sites().len(),
            frac.len()
        )));
    }
    let sites = structure
        .sites()
        .iter()
        .zip(frac)
        .map(|(site, point)| site.clone().with_frac(*point))
        .collect::<Result<Vec<_>, _>>()?;
    rebuild(structure, lattice, sites)
}
