//! Symmetry detection behind our own dataset type (Task 2.1.3).

use moyo::MoyoDataset;
use moyo::base::AngleTolerance;
use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure, SymOp};

use crate::convert::to_moyo;
use crate::spacegroup::setting;

/// Which standardized setting the search should prefer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SettingChoice {
    /// The setting spglib chooses (smallest Hall number).
    #[default]
    Spglib,
    /// Unique axis b, cell choice 1, hexagonal axes, origin choice 2.
    Standard,
    /// A specific Hall number from 1 to 530.
    HallNumber(u32),
}

/// Options for [`detect_symmetry`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SymmetryOptions {
    /// Position tolerance in Angstrom.
    pub symprec: f64,
    /// Angle tolerance in radians; `None` uses moyo's symprec-driven default.
    pub angle_tolerance: Option<f64>,
    /// Preferred standardized setting.
    pub setting: SettingChoice,
}

impl Default for SymmetryOptions {
    fn default() -> Self {
        Self {
            symprec: 1e-5,
            angle_tolerance: None,
            setting: SettingChoice::Spglib,
        }
    }
}

/// Result of a symmetry search.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SymmetryDataset {
    /// ITA space-group number.
    pub number: u32,
    /// Hall number of the detected setting.
    pub hall_number: u32,
    /// Hermann-Mauguin short symbol.
    pub hm_symbol: String,
    /// Pearson symbol of the standardized cell.
    pub pearson_symbol: String,
    /// Operations in the input cell.
    pub operations: Vec<SymOp>,
    /// Representative atom index for each input atom.
    pub orbits: Vec<usize>,
    /// Wyckoff letter for each input atom.
    pub wyckoffs: Vec<char>,
    /// Site-symmetry symbol for each input atom.
    pub site_symmetry_symbols: Vec<String>,
    /// Standardized lattice rows.
    pub std_lattice: [[f64; 3]; 3],
    /// Standardized fractional positions.
    pub std_positions: Vec<[f64; 3]>,
    /// Element symbol of each standardized position.
    pub std_elements: Vec<String>,
    /// Linear part of the input-to-standard transformation.
    pub std_linear: [[f64; 3]; 3],
    /// Primitive standardized lattice rows.
    pub prim_std_lattice: [[f64; 3]; 3],
    /// Primitive standardized fractional positions.
    pub prim_std_positions: Vec<[f64; 3]>,
    /// Element symbol of each primitive standardized position.
    pub prim_std_elements: Vec<String>,
    /// Mapping from input atoms to primitive standardized positions.
    pub mapping_std_prim: Vec<usize>,
    /// Origin shift from the input cell to the standardized cell.
    pub std_origin_shift: [f64; 3],
    /// Proper Cartesian rotation applied after lattice refinement.
    pub std_rotation_matrix: [[f64; 3]; 3],
}

impl SymmetryDataset {
    /// Number of atoms in the input cell.
    #[must_use]
    pub fn n_atoms(&self) -> usize {
        self.orbits.len()
    }

    /// Wyckoff letter of one input atom.
    #[must_use]
    pub fn wyckoff(&self, index: usize) -> Option<char> {
        self.wyckoffs.get(index).copied()
    }

    /// Indices of the atoms equivalent to `index`, including itself.
    #[must_use]
    pub fn equivalents(&self, index: usize) -> Vec<usize> {
        let Some(representative) = self.orbits.get(index) else {
            return Vec::new();
        };
        self.orbits
            .iter()
            .enumerate()
            .filter(|(_, orbit)| *orbit == representative)
            .map(|(atom, _)| atom)
            .collect()
    }
}

/// Detects the symmetry of a structure.
///
/// The structure is expanded to P1 first; the input is not modified.
///
/// # Errors
///
/// Returns [`Error::Symmetry`] when the search fails, the tolerance is
/// invalid or the structure is empty.
pub fn detect_symmetry(
    structure: &Structure,
    options: SymmetryOptions,
) -> Result<SymmetryDataset, Error> {
    if !options.symprec.is_finite() || options.symprec <= 0.0 {
        return Err(Error::Symmetry(format!(
            "symprec must be positive and finite, got {}",
            options.symprec
        )));
    }
    if let Some(angle) = options.angle_tolerance
        && (!angle.is_finite() || angle <= 0.0)
    {
        return Err(Error::Symmetry(format!(
            "angle tolerance must be positive and finite, got {angle}"
        )));
    }
    let p1 = structure.p1()?;
    let moyo_cell = to_moyo(&p1)?;
    let angle = options
        .angle_tolerance
        .map_or(AngleTolerance::Default, AngleTolerance::Radian);
    let dataset = MoyoDataset::new(
        &moyo_cell.cell,
        options.symprec,
        angle,
        setting(options.setting),
        false,
    )
    .map_err(|error| {
        Error::Symmetry(format!(
            "moyo search failed at symprec={} A: {error}",
            options.symprec
        ))
    })?;
    let positions: Vec<[f64; 3]> = p1.atoms.iter().map(|atom| atom.frac).collect();
    dataset_from_moyo(&dataset, &moyo_cell, &positions, options.symprec)
}

fn dataset_from_moyo(
    dataset: &MoyoDataset,
    moyo_cell: &crate::convert::MoyoCell,
    positions: &[[f64; 3]],
    symprec: f64,
) -> Result<SymmetryDataset, Error> {
    let n_atoms = positions.len();
    let operations = crate::convert::operations_from_moyo(&dataset.operations)?;
    let orbits = lift_orbits(dataset, n_atoms, symprec)?;
    let std_lattice = dataset.std_cell.lattice.basis.transpose().into();
    let std_positions = dataset.std_cell.positions_as_arrays();
    let std_elements = elements_for(&dataset.std_cell, moyo_cell)?;
    let prim_std_positions = dataset.prim_std_cell.positions_as_arrays();
    let prim_std_elements = elements_for(&dataset.prim_std_cell, moyo_cell)?;
    let std_linear: [[f64; 3]; 3] = dataset.std_linear.into();
    let prim_std_lattice = dataset.prim_std_cell.lattice.basis.transpose().into();
    let mut std_origin_shift = [0.0; 3];
    for (slot, value) in std_origin_shift
        .iter_mut()
        .zip(dataset.std_origin_shift.iter())
    {
        *slot = *value;
    }
    let std_rotation_matrix: [[f64; 3]; 3] = dataset.std_rotation_matrix.into();
    let wyckoffs = crate::wyckoff_origin::input_origin_letters(
        &operations,
        dataset.hall_number as u32,
        std_linear,
        positions,
        &dataset.wyckoffs,
    );
    Ok(SymmetryDataset {
        number: dataset.number as u32,
        hall_number: dataset.hall_number as u32,
        hm_symbol: dataset.hm_symbol.clone(),
        pearson_symbol: dataset.pearson_symbol.clone(),
        operations,
        orbits,
        wyckoffs,
        site_symmetry_symbols: dataset.site_symmetry_symbols.clone(),
        std_lattice,
        std_positions,
        std_elements,
        std_linear,
        prim_std_lattice,
        prim_std_positions,
        prim_std_elements,
        mapping_std_prim: dataset.mapping_std_prim.clone(),
        std_origin_shift,
        std_rotation_matrix,
    })
}

/// moyo reports orbits and Wyckoff letters for the input cell when the input
/// basis is the conventional cell; the lengths are checked here and a mismatch
/// is a hard error rather than silently misaligned data.
fn lift_orbits(dataset: &MoyoDataset, n_atoms: usize, symprec: f64) -> Result<Vec<usize>, Error> {
    if dataset.orbits.len() == n_atoms {
        return Ok(dataset.orbits.clone());
    }
    Err(Error::Symmetry(format!(
        "moyo returned {} orbits for {n_atoms} atoms at symprec={symprec} A; \
         the input cell is not the conventional cell",
        dataset.orbits.len()
    )))
}

fn elements_for(
    cell: &moyo::base::Cell,
    moyo_cell: &crate::convert::MoyoCell,
) -> Result<Vec<String>, Error> {
    cell.numbers
        .iter()
        .map(|number| {
            let index = usize::try_from(*number)
                .ok()
                .and_then(|value| value.checked_sub(1))
                .ok_or_else(|| Error::Symmetry(format!("invalid moyo atom type {number}")))?;
            let atom_type = moyo_cell
                .type_of_number
                .get(index)
                .ok_or_else(|| Error::Symmetry(format!("unknown moyo atom type {number}")))?;
            Ok(xtalforge_core::elements::element_by_z(atom_type.z)?
                .symbol
                .clone())
        })
        .collect()
}
