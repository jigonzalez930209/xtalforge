//! Maximal subgroup enumeration through moyo (Task 2.7.1).

use moyo::data::{HallNumber, Setting};
use serde::Serialize;
use xtalforge_core::{Error, SymOp};

use crate::convert::operations_from_moyo;

/// Subgroup family.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum SubgroupKind {
    /// Same translation lattice, smaller point group.
    Translationengleiche,
    /// Same point group, finite-index translation sublattice.
    Klassengleiche,
}

/// One subgroup conjugacy class representative.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Subgroup {
    /// Family.
    pub kind: SubgroupKind,
    /// Hall number of the parent.
    pub parent_hall_number: u32,
    /// Subgroup index in the parent.
    pub index: usize,
    /// Shortest covering-chain depth (translationengleiche only).
    pub depth: usize,
    /// Operations in the parent primitive basis.
    pub operations: Vec<SymOp>,
    /// Number of conjugates in the class.
    pub conjugates: usize,
    /// Identified space-group number, when moyo could name the subgroup.
    pub subgroup_number: Option<u32>,
    /// Identified Hall number, when moyo could name the subgroup.
    pub subgroup_hall_number: Option<u32>,
}

/// Primitive operations of a Hall setting.
///
/// # Errors
///
/// Returns [`Error::Symmetry`] for an unknown Hall number or an invalid
/// operation.
pub fn primitive_operations(hall_number: u32) -> Result<Vec<SymOp>, Error> {
    if !(1..=530).contains(&hall_number) {
        return Err(Error::Symmetry(format!(
            "Hall number must be in 1..=530, got {hall_number}"
        )));
    }
    let symbol = moyo::data::HallSymbol::from_hall_number(hall_number as HallNumber)
        .ok_or_else(|| Error::Symmetry(format!("unknown Hall symbol for {hall_number}")))?;
    operations_from_moyo(&symbol.primitive_traverse())
}

/// Translationengleiche subgroups, including the parent and the pure
/// translation subgroup.
///
/// # Errors
///
/// Returns [`Error::Symmetry`] when moyo cannot enumerate the subgroups.
pub fn translationengleiche_subgroups(
    hall_number: u32,
    epsilon: f64,
) -> Result<Vec<Subgroup>, Error> {
    let operations = primitive_operations(hall_number)?;
    let moyo_operations = to_moyo_operations(&operations)?;
    let classes =
        moyo::subgroup::enumerate_translationengleiche_subgroups(&moyo_operations, epsilon)
            .map_err(|error| Error::Symmetry(format!("subgroup enumeration failed: {error}")))?;
    let mut groups = Vec::new();
    for class in classes {
        let representative = class.representative;
        let operations = operations_from_moyo(&representative.operations)?;
        let identification = identify_space_group(&operations).ok();
        groups.push(Subgroup {
            kind: SubgroupKind::Translationengleiche,
            parent_hall_number: hall_number,
            index: representative.translationengleiche_index,
            depth: representative.depth,
            operations,
            conjugates: class.conjugates.len(),
            subgroup_number: identification.map(|value| value.0),
            subgroup_hall_number: identification.map(|value| value.1),
        });
    }
    Ok(groups)
}

/// Klassengleiche subgroups for a translation sublattice.
///
/// `transformation` is a row-basis matrix in the parent primitive basis; its
/// determinant must be positive and preserved by every parent rotation.
///
/// # Errors
///
/// Returns [`Error::Symmetry`] when moyo cannot enumerate the subgroups.
pub fn klassengleiche_subgroups(
    hall_number: u32,
    transformation: [[i32; 3]; 3],
    epsilon: f64,
) -> Result<Vec<Subgroup>, Error> {
    let operations = primitive_operations(hall_number)?;
    let moyo_operations = to_moyo_operations(&operations)?;
    let linear = nalgebra::Matrix3::from_row_slice(&[
        transformation[0][0],
        transformation[0][1],
        transformation[0][2],
        transformation[1][0],
        transformation[1][1],
        transformation[1][2],
        transformation[2][0],
        transformation[2][1],
        transformation[2][2],
    ]);
    let classes =
        moyo::subgroup::enumerate_klassengleiche_subgroups(&moyo_operations, &linear, epsilon)
            .map_err(|error| Error::Symmetry(format!("subgroup enumeration failed: {error}")))?;
    let mut groups = Vec::new();
    for class in classes {
        let representative = class.representative;
        let operations = operations_from_moyo(&representative.parent_operations)?;
        let identification = identify_space_group(&operations).ok();
        groups.push(Subgroup {
            kind: SubgroupKind::Klassengleiche,
            parent_hall_number: hall_number,
            index: representative.klassengleiche_index,
            depth: 1,
            operations,
            conjugates: class.conjugates.len(),
            subgroup_number: identification.map(|value| value.0),
            subgroup_hall_number: identification.map(|value| value.1),
        });
    }
    Ok(groups)
}

/// Identifies the space group of a primitive operation list.
///
/// # Errors
///
/// Returns [`Error::Symmetry`] when moyo cannot identify the group.
pub fn identify_space_group(operations: &[SymOp]) -> Result<(u32, u32), Error> {
    let moyo_operations = to_moyo_operations(operations)?;
    let group = moyo::identify::SpaceGroup::new(&moyo_operations, Setting::Standard, 1e-5)
        .map_err(|error| Error::Symmetry(format!("space-group identification failed: {error}")))?;
    Ok((group.number as u32, group.hall_number as u32))
}

fn to_moyo_operations(operations: &[SymOp]) -> Result<Vec<moyo::base::Operation>, Error> {
    use nalgebra::{Matrix3, Vector3};
    operations
        .iter()
        .map(|operation| {
            let rotation = operation.rotation();
            let translation = operation.translation();
            Ok(moyo::base::Operation::new(
                Matrix3::from_row_slice(&[
                    rotation[0][0],
                    rotation[0][1],
                    rotation[0][2],
                    rotation[1][0],
                    rotation[1][1],
                    rotation[1][2],
                    rotation[2][0],
                    rotation[2][1],
                    rotation[2][2],
                ]),
                Vector3::new(translation[0], translation[1], translation[2]),
            ))
        })
        .collect()
}
