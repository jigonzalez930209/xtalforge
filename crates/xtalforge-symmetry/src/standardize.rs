//! Standardized and primitive cells and idealized positions (Task 2.5).

use xtalforge_core::{Error, Site, Structure, StructureBuilder};

use crate::dataset::{SymmetryDataset, SymmetryOptions, detect_symmetry};
use crate::spacegroup::SpaceGroup;
use crate::workflow::detect_and_symmetrize;

/// A standardized structure and the transformation that produced it.
#[derive(Debug, Clone)]
pub struct StandardizedStructure {
    /// Asymmetric-unit structure in the standardized setting.
    pub structure: Structure,
    /// Origin shift from the input cell.
    pub origin_shift: [f64; 3],
    /// Proper Cartesian rotation applied after lattice refinement.
    pub rotation_matrix: [[f64; 3]; 3],
    /// The standardized dataset.
    pub dataset: SymmetryDataset,
}

/// Standardizes a structure: conventional cell, detected setting, refined
/// positions and the asymmetric unit.
///
/// # Errors
///
/// Returns the errors of the search or a builder error.
pub fn standardize(
    structure: &Structure,
    options: SymmetryOptions,
) -> Result<StandardizedStructure, Error> {
    let dataset = detect_symmetry(structure, options)?;
    let group = SpaceGroup::from_hall_number(dataset.hall_number)?;
    let temp = full_cell_structure(
        structure,
        dataset.std_lattice,
        &dataset.std_positions,
        &dataset.std_elements,
    )?;
    let reduced = detect_and_symmetrize(
        &temp,
        SymmetryOptions {
            symprec: options.symprec,
            ..SymmetryOptions::default()
        },
    )?;
    let mut std_structure = reduced.structure;
    // Keep the requested Hall setting and the input metadata.
    if std_structure.hall() != group.hall_symbol {
        std_structure = std_structure_with_hall(std_structure, &group)?;
    }
    Ok(StandardizedStructure {
        structure: std_structure,
        origin_shift: dataset.std_origin_shift,
        rotation_matrix: dataset.std_rotation_matrix,
        dataset,
    })
}

/// Builds the standardized primitive cell and its asymmetric unit.
///
/// # Errors
///
/// Returns the errors of the search or a builder error.
pub fn primitive_cell(
    structure: &Structure,
    options: SymmetryOptions,
) -> Result<StandardizedStructure, Error> {
    let dataset = detect_symmetry(structure, options)?;
    let temp = full_cell_structure(
        structure,
        dataset.prim_std_lattice,
        &dataset.prim_std_positions,
        &dataset.prim_std_elements,
    )?;
    let reduced = detect_and_symmetrize(
        &temp,
        SymmetryOptions {
            symprec: options.symprec,
            ..SymmetryOptions::default()
        },
    )?;
    Ok(StandardizedStructure {
        structure: reduced.structure,
        origin_shift: dataset.std_origin_shift,
        rotation_matrix: dataset.std_rotation_matrix,
        dataset,
    })
}

/// Idealizes positions by returning the standardized, refined structure.
///
/// moyo refines positions to the detected symmetry, so this is the same
/// operation as [`standardize`].
///
/// # Errors
///
/// Returns the errors of [`standardize`].
pub fn idealize(
    structure: &Structure,
    options: SymmetryOptions,
) -> Result<StandardizedStructure, Error> {
    standardize(structure, options)
}

fn full_cell_structure(
    source: &Structure,
    lattice: [[f64; 3]; 3],
    positions: &[[f64; 3]],
    elements: &[String],
) -> Result<Structure, Error> {
    let lattice = xtalforge_core::Lattice::new(lattice)?;
    let mut builder = StructureBuilder::new(source.name(), lattice);
    for (index, (position, element)) in positions.iter().zip(elements.iter()).enumerate() {
        // The standardized cell can duplicate a site; preserve its typed
        // properties when the element is unambiguous. A multi-site element
        // keeps the geometry but falls back to a fresh site because moyo does
        // not expose a site-to-site mapping for the standardized cell.
        let matches: Vec<&Site> = source
            .sites()
            .iter()
            .filter(|site| site.element() == element)
            .collect();
        let site = if matches.len() == 1 {
            let Some(site) = matches.first() else {
                continue;
            };
            (*site)
                .clone()
                .with_label(format!("{element}{}", index + 1))?
                .with_frac(*position)?
        } else {
            Site::new(format!("{element}{}", index + 1), element, *position)?
        };
        builder = builder.site(site);
    }
    builder.build()
}

fn std_structure_with_hall(structure: Structure, group: &SpaceGroup) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new(structure.name(), *structure.lattice())
        .hall(group.hall_symbol.clone());
    for operation in &group.operations {
        builder = builder.operation(*operation);
    }
    for site in structure.sites() {
        builder = builder.site(site.clone());
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    builder.build()
}

/// Changes the setting of a structure to a specific Hall number.
///
/// The structure is standardized in the requested setting; the physical
/// geometry is preserved while the origin, axes and cell choice change.
///
/// # Errors
///
/// Returns [`Error::Symmetry`] for an unknown Hall number or the errors of
/// [`standardize`].
pub fn change_setting(
    structure: &Structure,
    hall_number: u32,
) -> Result<StandardizedStructure, Error> {
    if !(1..=530).contains(&hall_number) {
        return Err(Error::Symmetry(format!(
            "Hall number must be in 1..=530, got {hall_number}"
        )));
    }
    standardize(
        structure,
        SymmetryOptions {
            setting: crate::SettingChoice::HallNumber(hall_number),
            ..SymmetryOptions::default()
        },
    )
}
