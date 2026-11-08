//! Guest placement (Task 6.4).

use xtalforge_core::{Error, Site, Structure, StructureBuilder};

/// Places a guest molecule at a Cartesian offset inside a host cell.
///
/// Every guest atom must keep at least `clearance` Angstrom from the host
/// atoms (periodic minimum image). The guest coordinates are relative to its
/// own origin.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty guest, a bad clearance or a
/// clash.
pub fn place_guest(
    host: &Structure,
    guest: &[(String, [f64; 3])],
    offset: [f64; 3],
    clearance: f64,
) -> Result<Structure, Error> {
    if guest.is_empty() {
        return Err(Error::InvalidInput(
            "a guest needs at least one atom".to_owned(),
        ));
    }
    if !clearance.is_finite() || clearance <= 0.0 {
        return Err(Error::InvalidInput(
            "guest clearance must be positive".to_owned(),
        ));
    }
    let p1 = host.p1()?;
    let mut builder =
        StructureBuilder::new(host.name(), *host.lattice()).hall(host.hall().to_owned());
    for operation in host.operations() {
        builder = builder.operation(*operation);
    }
    for site in host.sites() {
        builder = builder.site(site.clone());
    }
    for (index, (element, point)) in guest.iter().enumerate() {
        let cart = [
            point[0] + offset[0],
            point[1] + offset[1],
            point[2] + offset[2],
        ];
        let frac = host.lattice().to_frac(cart);
        let nearest = p1
            .atoms
            .iter()
            .map(|atom| p1.lattice.distance(frac, atom.frac))
            .fold(f64::INFINITY, f64::min);
        if nearest < clearance {
            return Err(Error::InvalidInput(format!(
                "guest atom {index} clashes with the host at {nearest:.3} A"
            )));
        }
        builder = builder.site(Site::new(
            format!("{element}g{}", index + 1),
            element,
            frac,
        )?);
    }
    for (key, value) in host.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}
