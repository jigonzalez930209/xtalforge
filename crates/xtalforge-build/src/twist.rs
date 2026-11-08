//! Commensurate graphene twists (Task 6.5.3).

use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

use crate::carbon::graphene;

/// Builds two graphene layers with a relative twist.
///
/// Only 0 and 60 degrees are exact in the small cell built here; any other
/// angle is approximated by the nearest of the two and the choice is reported
/// in the structure metadata (`twist_approximation`).
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for zero repeats/layers or a non-finite
/// angle, or the errors of [`graphene`].
pub fn commensurate_twist(
    n: u32,
    m: u32,
    theta_degrees: f64,
    layers: u32,
) -> Result<Structure, Error> {
    if n == 0 || m == 0 || layers == 0 || !theta_degrees.is_finite() {
        return Err(Error::InvalidInput(
            "twist repeats, layers and angle must be valid".to_owned(),
        ));
    }
    let nearest = if theta_degrees.rem_euclid(60.0) <= 30.0 {
        0.0
    } else {
        60.0
    };
    let base = graphene(n, m, layers, 3.35, 15.0)?;
    let matrix = base.lattice().matrix();
    let lattice = Lattice::new(matrix)?;
    let mut builder = StructureBuilder::new("twisted-graphene", lattice)
        .metadata("twist_requested", serde_json::json!(theta_degrees))
        .metadata("twist_approximation", serde_json::json!(nearest));
    for (index, site) in base.sites().iter().enumerate() {
        let frac = site.frac();
        let shifted = if nearest > 0.0 {
            [
                (frac[0] + 1.0 / 3.0).rem_euclid(1.0),
                (frac[1] + 1.0 / 3.0).rem_euclid(1.0),
                frac[2],
            ]
        } else {
            frac
        };
        builder = builder.site(Site::new(
            format!("C{}", index + 1),
            site.element(),
            shifted,
        )?);
    }
    builder.build()
}
