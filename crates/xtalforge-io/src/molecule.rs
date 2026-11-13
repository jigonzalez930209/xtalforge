//! Helpers shared by the molecular file formats (Task 13.2.5).
//!
//! A molecule has no cell in MOL2 or SDF, so the readers place it in a padded
//! P1 box exactly like the plain-XYZ path. Labels and explicit bonds are kept.

use xtalforge_core::{Error, ExplicitBond, Lattice, Site, Structure, StructureBuilder};

/// Padding around a molecule without a cell, in Angstrom.
pub(crate) const PADDING: f64 = 5.0;

/// One atom read from a molecular file.
pub(crate) struct MoleculeAtom {
    /// Site label.
    pub label: String,
    /// Element symbol or atom type.
    pub element: String,
    /// Cartesian coordinates in Angstrom.
    pub cart: [f64; 3],
}

/// Builds a P1 molecule in a padded orthorhombic box.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for no atoms, a non-finite coordinate or an
/// unknown element.
pub(crate) fn molecule_box(
    name: &str,
    atoms: &[MoleculeAtom],
    bonds: &[ExplicitBond],
) -> Result<Structure, Error> {
    if atoms.is_empty() {
        return Err(Error::InvalidInput(
            "a molecule file needs at least one atom".to_owned(),
        ));
    }
    let mut minimum = [f64::INFINITY; 3];
    let mut maximum = [f64::NEG_INFINITY; 3];
    for atom in atoms {
        if !atom.cart.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "molecule coordinates must be finite".to_owned(),
            ));
        }
        for (axis, value) in atom.cart.iter().enumerate() {
            if let (Some(low), Some(high)) = (minimum.get_mut(axis), maximum.get_mut(axis)) {
                *low = low.min(*value);
                *high = high.max(*value);
            }
        }
    }
    let lengths = [
        maximum[0] - minimum[0] + 2.0 * PADDING,
        maximum[1] - minimum[1] + 2.0 * PADDING,
        maximum[2] - minimum[2] + 2.0 * PADDING,
    ];
    let lattice = Lattice::orthorhombic(lengths[0], lengths[1], lengths[2])?;
    let origin = [
        minimum[0] - PADDING,
        minimum[1] - PADDING,
        minimum[2] - PADDING,
    ];
    let mut builder = StructureBuilder::new(name, lattice);
    for atom in atoms {
        let frac = [
            (atom.cart[0] - origin[0]) / lengths[0],
            (atom.cart[1] - origin[1]) / lengths[1],
            (atom.cart[2] - origin[2]) / lengths[2],
        ];
        builder = builder.site(Site::new(atom.label.clone(), &atom.element, frac)?);
    }
    for bond in bonds {
        if bond.i < atoms.len() && bond.j < atoms.len() {
            builder = builder.explicit_bond(*bond);
        }
    }
    builder.build()
}

/// Parses a numeric field, reporting the format in the error.
pub(crate) fn parse_number(text: &str, what: &str) -> Result<f64, Error> {
    text.trim()
        .parse()
        .map_err(|_| Error::InvalidInput(format!("bad {what} number: {text:?}")))
}

/// Parses an integer field, reporting the format in the error.
pub(crate) fn parse_index(text: &str, what: &str) -> Result<usize, Error> {
    text.trim()
        .parse()
        .map_err(|_| Error::InvalidInput(format!("bad {what} index: {text:?}")))
}

/// Makes a label unique by appending `_2`, `_3`, ... when it repeats.
pub(crate) fn unique_label(label: &str, used: &mut std::collections::BTreeSet<String>) -> String {
    if used.insert(label.to_owned()) {
        return label.to_owned();
    }
    let mut suffix = 2usize;
    loop {
        let candidate = format!("{label}_{suffix}");
        if used.insert(candidate.clone()) {
            return candidate;
        }
        suffix += 1;
    }
}
