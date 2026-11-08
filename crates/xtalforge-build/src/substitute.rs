//! Clearance checks and site substitution (Task 6.3).

use xtalforge_core::{Error, Site, Structure, StructureBuilder};

/// Minimum distance from a fractional point to any atom centre.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-finite point.
pub fn clearance(structure: &Structure, point: [f64; 3]) -> Result<f64, Error> {
    if !point.iter().all(|value| value.is_finite()) {
        return Err(Error::InvalidInput(
            "clearance point must be finite".to_owned(),
        ));
    }
    let p1 = structure.p1()?;
    Ok(p1
        .atoms
        .iter()
        .map(|atom| p1.lattice.distance(point, atom.frac))
        .fold(f64::INFINITY, f64::min))
}

/// Replaces the element of one asymmetric site.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an out-of-range index, an unknown
/// element or a builder failure.
pub fn substitute(
    structure: &Structure,
    site_index: usize,
    element: &str,
) -> Result<Structure, Error> {
    let site = structure.sites().get(site_index).ok_or_else(|| {
        Error::InvalidInput(format!(
            "site {site_index} is outside the structure ({} sites)",
            structure.sites().len()
        ))
    })?;
    let mut builder = StructureBuilder::new(structure.name(), *structure.lattice())
        .hall(structure.hall().to_owned());
    for operation in structure.operations() {
        builder = builder.operation(*operation);
    }
    for (index, current) in structure.sites().iter().enumerate() {
        if index == site_index {
            let replacement = Site::new(current.label(), element, current.frac())?
                .with_occupancy(current.occupancy())?;
            builder = builder.site(replacement);
        } else {
            builder = builder.site(current.clone());
        }
    }
    let _ = site;
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}
