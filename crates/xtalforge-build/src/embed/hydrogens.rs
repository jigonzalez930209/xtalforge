//! Hydrogen inference from SMILES bond orders.
//!
//! The parser records aromatic bonds as order `1.5`, so the graph carries
//! enough information to rebuild the hydrogens a SMILES string leaves
//! implicit. Each covered atom receives the difference between its target
//! valence and its rounded bond-order sum, placed in the directions the
//! already embedded heavy atoms leave free.

use xtalforge_core::Error;

use crate::hydrogens::geometry::directions_from_points;
use crate::hydrogens::{inference_geometry, inference_missing};
use crate::smiles::SmilesMolecule;

/// Hydrogens implied by the graph, in atom order, as element and point.
///
/// The heavy atoms keep their graph order; every inferred hydrogen follows
/// them, grouped by the atom that owns it.
pub(super) fn inferred_hydrogens(
    graph: &SmilesMolecule,
    points: &[[f64; 3]],
) -> Result<Vec<(String, [f64; 3])>, Error> {
    let mut hydrogens = Vec::new();
    for (index, atom) in graph.atoms.iter().enumerate() {
        let Some(origin) = points.get(index).copied() else {
            return Err(Error::InvalidInput(format!(
                "embedding point {index} is missing"
            )));
        };
        let mut neighbours = Vec::new();
        let mut used_order = 0.0;
        for bond in &graph.bonds {
            let other = if bond.a == index {
                bond.b
            } else if bond.b == index {
                bond.a
            } else {
                continue;
            };
            used_order += bond.order;
            let Some(point) = points.get(other) else {
                return Err(Error::InvalidInput(format!(
                    "embedding point {other} is missing"
                )));
            };
            neighbours.push(*point);
        }
        let Some(missing) = inference_missing(&atom.element, atom.charge, used_order) else {
            continue;
        };
        if missing == 0 {
            continue;
        }
        let (length, angle) = inference_geometry(&atom.element);
        for direction in directions_from_points(origin, &neighbours, missing, angle)? {
            hydrogens.push((
                "H".to_owned(),
                [
                    origin[0] + direction[0] * length,
                    origin[1] + direction[1] * length,
                    origin[2] + direction[2] * length,
                ],
            ));
        }
    }
    Ok(hydrogens)
}
