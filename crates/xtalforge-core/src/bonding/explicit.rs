//! Explicit bonds expanded over every symmetry image (Task 14.3.2).
//!
//! A stored bond names two atoms of the P1 cell and the lattice image of the
//! second. Drawing one Fe-O bond must draw the whole orbit, exactly as VESTA
//! does, so every operation of the space group is applied to both ends at
//! once. Applying the operation to both ends is what keeps a site on a special
//! position from producing one bond where the symmetry demands several.

use crate::{ExplicitBond, P1Cell, Structure, SymOp};

use super::canonical;

/// Tolerance in Angstrom for finding the P1 atom that an operation generates.
pub(crate) const IMAGE_TOLERANCE: f64 = 1e-3;

/// Expands one stored bond into every symmetry image it draws.
///
/// The input indices name atoms of the P1 cell. Every operation of the
/// structure is applied to both ends; the generated points are matched to the
/// P1 atoms that are stored there and the bond image is read back from the
/// matched atoms, so a partner that wrapped across the cell face keeps the
/// translation that actually joins the two atoms. Duplicate spellings of the
/// same edge are removed and the result is ordered by `(i, j, image)`.
///
/// A structure with no operations expands with the identity only.
#[must_use]
pub fn explicit_bond_images(
    structure: &Structure,
    p1: &P1Cell,
    bond: &ExplicitBond,
) -> Vec<ExplicitBond> {
    // The CIF reader and the writer store explicit bonds with P1 indices
    // (`_xtal_bond_atom_1`) and mark the structure; the editor stores them on
    // the asymmetric unit. Expanding a P1 bond again would paint every
    // symmetry image onto the same contact.
    if structure
        .metadata()
        .get("explicit_bonds_p1")
        .and_then(serde_json::Value::as_bool)
        .unwrap_or(false)
    {
        return vec![*bond];
    }
    let identity = [SymOp::identity()];
    let operations = if structure.operations().is_empty() {
        identity.as_slice()
    } else {
        structure.operations()
    };
    let Some(first) = p1.atoms.get(bond.i) else {
        return Vec::new();
    };
    let Some(second) = p1.atoms.get(bond.j) else {
        return Vec::new();
    };
    let far = [
        second.frac[0] + f64::from(bond.image[0]),
        second.frac[1] + f64::from(bond.image[1]),
        second.frac[2] + f64::from(bond.image[2]),
    ];
    let mut found: Vec<ExplicitBond> = Vec::new();
    for operation in operations {
        let here_raw = operation.apply(first.frac);
        let (here, _) = wrap_with_tau(here_raw);
        let there_raw = operation.apply(far);
        let (there, _) = wrap_with_tau(there_raw);
        let (Some(a), Some(b)) = (nearest_atom(p1, here), nearest_atom(p1, there)) else {
            continue;
        };
        let (Some(atom_a), Some(atom_b)) = (p1.atoms.get(a), p1.atoms.get(b)) else {
            continue;
        };
        // Read the image from the atoms that were found, not from the wrapped
        // points: an atom stored at 0 is found for a point wrapped to 0.99999,
        // and taking the image from the point would hang the bond on the copy a
        // cell away.
        let image = [
            round(there_raw[0] - atom_b.frac[0] - (here_raw[0] - atom_a.frac[0])),
            round(there_raw[1] - atom_b.frac[1] - (here_raw[1] - atom_a.frac[1])),
            round(there_raw[2] - atom_b.frac[2] - (here_raw[2] - atom_a.frac[2])),
        ];
        let (i, j, image) = canonical(a, b, image);
        let expanded = ExplicitBond {
            i,
            j,
            image,
            ..*bond
        };
        if !found
            .iter()
            .any(|existing| same_record(existing, &expanded))
        {
            found.push(expanded);
        }
    }
    found.sort_by_key(|bond| (bond.i, bond.j, bond.image, bond.remove));
    found
}

fn same_record(first: &ExplicitBond, second: &ExplicitBond) -> bool {
    first.i == second.i
        && first.j == second.j
        && first.image == second.image
        && first.remove == second.remove
}

fn nearest_atom(p1: &P1Cell, frac: [f64; 3]) -> Option<usize> {
    let mut best: Option<(f64, usize)> = None;
    for (index, atom) in p1.atoms.iter().enumerate() {
        let distance = p1.lattice.distance(atom.frac, frac);
        if distance <= IMAGE_TOLERANCE && best.is_none_or(|(value, _)| distance < value) {
            best = Some((distance, index));
        }
    }
    best.map(|(_, index)| index)
}

fn wrap_with_tau(raw: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    let floor = raw.map(f64::floor);
    (
        [raw[0] - floor[0], raw[1] - floor[1], raw[2] - floor[2]],
        [-floor[0], -floor[1], -floor[2]],
    )
}

fn round(value: f64) -> i32 {
    value.round() as i32
}
