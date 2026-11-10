//! Detection of a file that lists the whole cell under a non-P1 group.

use xtalforge_core::{Structure, neighbors};

/// Distance in Angstrom under which two images are the same position (the P1
/// expansion tolerance).
const COINCIDENT_TOLERANCE: f64 = 0.05;

/// A warning when atoms of different sites coincide in the P1 cell and the
/// occupancies add up to more than one.
///
/// A file that lists every atom of the cell and also names a space group makes
/// each orbit appear several times; real disorder components have partial
/// occupancy, so their total stays at one.
pub(super) fn warning(structure: &Structure) -> Option<String> {
    if structure.is_p1() {
        return None;
    }
    let p1 = structure.p1().ok()?;
    let pairs = neighbors(structure, COINCIDENT_TOLERANCE).ok()?;
    let count = pairs
        .iter()
        .filter(|pair| pair.image == [0, 0, 0])
        .filter(|pair| {
            let (Some(first), Some(second)) = (p1.atoms.get(pair.i), p1.atoms.get(pair.j)) else {
                return false;
            };
            first.site_idx != second.site_idx
                && first.element == second.element
                && first.occupancy + second.occupancy > 1.0 + 1e-6
        })
        .count();
    (count > 0).then(|| {
        format!(
            "{count} pairs of atoms from different sites coincide in the P1 cell with total occupancy above one; the file may list the whole cell under a non-P1 group"
        )
    })
}
