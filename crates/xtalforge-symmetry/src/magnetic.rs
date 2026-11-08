//! Collinear magnetic space groups (Task 2.8.2, experimental).

use moyo::MoyoMagneticDataset;
use moyo::base::{AngleTolerance, Collinear, MagneticCell, RotationMagneticMomentAction};
use serde::Serialize;
use xtalforge_core::{Error, Structure};

use crate::convert::to_moyo;

/// One magnetic operation.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct MagneticOperation {
    /// Integer rotation.
    pub rotation: [[i32; 3]; 3],
    /// Fractional translation.
    pub translation: [f64; 3],
    /// Whether time reversal is applied.
    pub time_reversal: bool,
}

/// Collinear magnetic symmetry dataset.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct MagneticDataset {
    /// UNI number of the magnetic space group (1 to 1651).
    pub uni_number: u32,
    /// Magnetic operations in the input cell.
    pub operations: Vec<MagneticOperation>,
    /// Representative atom index for each input atom.
    pub orbits: Vec<usize>,
}

/// Detects the collinear magnetic space group of a structure.
///
/// Magnetic moments come from the site property `magmom` (a number); a site
/// without it is treated as non-magnetic (`0`). The action on the moments is
/// axial.
///
/// # Errors
///
/// Returns [`Error::Symmetry`] when the search fails or the tolerance is
/// invalid.
pub fn detect_magnetic_symmetry(
    structure: &Structure,
    symprec: f64,
) -> Result<MagneticDataset, Error> {
    if !symprec.is_finite() || symprec <= 0.0 {
        return Err(Error::Symmetry(format!(
            "symprec must be positive and finite, got {symprec}"
        )));
    }
    let p1 = structure.p1()?;
    let moyo_cell = to_moyo(&p1)?;
    let moments: Vec<Collinear> = p1
        .atoms
        .iter()
        .map(|atom| {
            let moment = structure
                .sites()
                .get(atom.site_idx)
                .and_then(|site| site.props().get("magmom"))
                .and_then(serde_json::Value::as_f64)
                .unwrap_or(0.0);
            Collinear(moment)
        })
        .collect();
    let magnetic_cell = MagneticCell::from_cell(moyo_cell.cell, moments);
    let dataset = MoyoMagneticDataset::new(
        &magnetic_cell,
        symprec,
        AngleTolerance::Default,
        None,
        RotationMagneticMomentAction::Axial,
        false,
    )
    .map_err(|error| {
        Error::Symmetry(format!(
            "magnetic search failed at symprec={symprec} A: {error}"
        ))
    })?;
    let operations = dataset
        .magnetic_operations
        .iter()
        .map(|operation| MagneticOperation {
            rotation: operation.operation.rotation_as_array(),
            translation: operation.operation.translation_as_array(),
            time_reversal: operation.time_reversal,
        })
        .collect();
    Ok(MagneticDataset {
        uni_number: dataset.uni_number as u32,
        operations,
        orbits: dataset.orbits,
    })
}
