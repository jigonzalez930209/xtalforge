//! Conversions between the core structure model and moyo cells (Task 2.1.1).
//!
//! Atom-type policy: a type is the pair `(element, occupancy)`, with the
//! occupancy rounded to `OCCUPANCY_DECIMALS` decimals. A mixed Fe/Ni site is
//! therefore never equivalent to a pure Fe site, and a dummy `X` gets its own
//! type. moyo only accepts integer species numbers, so the pairs are mapped to
//! dense indices in first-appearance order.

use std::collections::BTreeMap;

use moyo::base::{Cell, Lattice, Position};
use xtalforge_core::elements::element_by_symbol;
use xtalforge_core::{Error, P1Cell, SymOp};

/// Decimal places used in the occupancy part of an atom type.
pub const OCCUPANCY_DECIMALS: i32 = 3;

/// Dense atom type used by the moyo conversion.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct AtomType {
    /// Atomic number, `0` for the dummy species.
    pub z: u8,
    /// Occupancy multiplied by `10^OCCUPANCY_DECIMALS` and rounded.
    pub occupancy: i64,
}

impl AtomType {
    /// Builds a type from an element symbol and occupancy.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] for an unknown element or a non-finite
    /// occupancy.
    pub fn new(symbol: &str, occupancy: f64) -> Result<Self, Error> {
        if !occupancy.is_finite() {
            return Err(Error::InvalidInput(
                "occupancy must be finite to build an atom type".to_owned(),
            ));
        }
        let z = element_by_symbol(symbol)?.z;
        let factor = 10f64.powi(OCCUPANCY_DECIMALS);
        Ok(Self {
            z,
            occupancy: (occupancy * factor).round() as i64,
        })
    }
}

/// A P1 cell ready for a moyo search.
pub(crate) struct MoyoCell {
    pub(crate) cell: Cell,
}

/// Converts a P1 cell to a moyo cell with dense atom types.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty cell, an unknown element or a
/// non-finite coordinate.
pub(crate) fn to_moyo(cell: &P1Cell) -> Result<MoyoCell, Error> {
    if cell.atoms.is_empty() {
        return Err(Error::Symmetry(
            "a symmetry search needs at least one atom".to_owned(),
        ));
    }
    let lattice = Lattice::from_basis(cell.lattice.matrix());
    let mut mapping: BTreeMap<AtomType, i32> = BTreeMap::new();
    let mut numbers = Vec::with_capacity(cell.atoms.len());
    let mut positions = Vec::with_capacity(cell.atoms.len());
    for atom in &cell.atoms {
        if !atom.frac.iter().all(|value| value.is_finite()) {
            return Err(Error::InvalidInput(
                "an atom has non-finite coordinates".to_owned(),
            ));
        }
        let atom_type = AtomType::new(&atom.element, atom.occupancy)?;
        let next = mapping.len() as i32 + 1;
        let number = *mapping.entry(atom_type.clone()).or_insert(next);
        numbers.push(number);
        positions.push(Position::new(atom.frac[0], atom.frac[1], atom.frac[2]));
    }
    Ok(MoyoCell {
        cell: Cell::new(lattice, positions, numbers),
    })
}

/// Converts a moyo operation to the core [`SymOp`].
///
/// # Errors
///
/// Returns [`Error::Symmetry`] if the rotation or translation cannot be
/// represented by the core type.
pub(crate) fn operation_from_moyo(operation: &moyo::base::Operation) -> Result<SymOp, Error> {
    SymOp::new(
        operation.rotation_as_array(),
        operation.translation_as_array(),
    )
    .map_err(|error| Error::Symmetry(format!("moyo returned an invalid operation: {error}")))
}

/// Converts the operation list of a moyo dataset.
///
/// # Errors
///
/// Returns the errors of [`operation_from_moyo`].
pub(crate) fn operations_from_moyo(
    operations: &[moyo::base::Operation],
) -> Result<Vec<SymOp>, Error> {
    operations
        .iter()
        .map(operation_from_moyo)
        .collect::<Result<Vec<_>, _>>()
}
