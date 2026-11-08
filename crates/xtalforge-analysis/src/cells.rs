//! Free cell parameters and Bravais extinction symbols (Tasks 7.5.3/7.5.4).
//!
//! A refinement form must offer only the numbers a space group leaves free:
//! offering six boxes for a tetragonal cell lets a fit start from
//! `a = 4.59, b = 4.60`, and the lattice ties then silently overwrite the
//! person's input. [`free_cell_parameters`] reports the numbers a caller may
//! refine for a space-group number in the standard International Tables
//! setting.
//!
//! # Trigonal settings
//!
//! Rhombohedral groups (`R3`, `R-3`, `R32`, `R3m`, `R3c`, `R-3m`, `R-3c`, and
//! the pairs below them) are reported in **rhombohedral axes** as `a` and
//! `alpha`, matching the reference `powder/cell.py`. Their hexagonal-axis
//! description is `a` and `c`; use [`free_cell_parameters_with_setting`] to ask
//! for it explicitly. P-lattice trigonal groups always report `a` and `c`.
//!
//! # Extinction symbols
//!
//! [`extinction_symbol`] returns the centring letter of the standard setting:
//! one of `P`, `A`, `B`, `C`, `I`, `F`, `R`. The 14 Bravais lattices use seven
//! letters because `A`, `B` and `C` are the same lattice in different
//! orientations; standard orthorhombic settings use `A` (groups 38-41) and `C`,
//! while `B` is accepted for non-standard settings but is never returned by the
//! standard-setting table.

use serde::{Deserialize, Serialize};
use xtalforge_core::Error;

/// The six cell numbers, in the order a cell is written.
pub const CELL_PARAMETER_NAMES: [&str; 6] = ["a", "b", "c", "alpha", "beta", "gamma"];

/// Axis choice for a trigonal cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum TrigonalSetting {
    /// Hexagonal axes: `a` and `c` free, `gamma = 120` degrees.
    HexagonalAxes,
    /// Rhombohedral axes: `a` and `alpha` free, all three angles equal.
    RhombohedralAxes,
}

/// Returns the cell numbers a space group leaves free.
///
/// The mapping follows the crystal system of the standard setting:
///
/// | Number | System | Free parameters |
/// |---|---|---|
/// | 1-2 | triclinic | `a b c alpha beta gamma` |
/// | 3-15 | monoclinic, unique axis b | `a b c beta` |
/// | 16-74 | orthorhombic | `a b c` |
/// | 75-142 | tetragonal | `a c` |
/// | 143-167 | trigonal P | `a c` |
/// | 143-167 | trigonal R | `a alpha` (rhombohedral axes) |
/// | 168-194 | hexagonal | `a c` |
/// | 195-230 | cubic | `a` |
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the number is outside `1..=230`.
pub fn free_cell_parameters(space_group_number: u32) -> Result<Vec<&'static str>, Error> {
    let setting = if is_rhombohedral(space_group_number) {
        TrigonalSetting::RhombohedralAxes
    } else {
        TrigonalSetting::HexagonalAxes
    };
    free_cell_parameters_with_setting(space_group_number, setting)
}

/// Returns the free cell numbers for an explicit trigonal setting.
///
/// For every crystal system other than trigonal the setting is ignored.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the number is outside `1..=230`, or
/// when [`TrigonalSetting::RhombohedralAxes`] is requested for a trigonal
/// P-lattice group (which has no rhombohedral description).
pub fn free_cell_parameters_with_setting(
    space_group_number: u32,
    setting: TrigonalSetting,
) -> Result<Vec<&'static str>, Error> {
    if !(1..=230).contains(&space_group_number) {
        return Err(Error::InvalidInput(format!(
            "space group number must be in 1..=230, got {space_group_number}"
        )));
    }
    let parameters = match space_group_number {
        1..=2 => vec!["a", "b", "c", "alpha", "beta", "gamma"],
        3..=15 => vec!["a", "b", "c", "beta"],
        16..=74 => vec!["a", "b", "c"],
        75..=142 => vec!["a", "c"],
        143..=167 => match setting {
            TrigonalSetting::HexagonalAxes => vec!["a", "c"],
            TrigonalSetting::RhombohedralAxes if is_rhombohedral(space_group_number) => {
                vec!["a", "alpha"]
            }
            TrigonalSetting::RhombohedralAxes => {
                return Err(Error::InvalidInput(format!(
                    "space group {space_group_number} is a trigonal P lattice and has no rhombohedral setting"
                )));
            }
        },
        168..=194 => vec!["a", "c"],
        _ => vec!["a"],
    };
    Ok(parameters)
}

/// Centring letter of the standard setting for a space-group number.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the number is outside `1..=230`.
pub fn extinction_symbol(space_group_number: u32) -> Result<&'static str, Error> {
    if !(1..=230).contains(&space_group_number) {
        return Err(Error::InvalidInput(format!(
            "space group number must be in 1..=230, got {space_group_number}"
        )));
    }
    Ok(match space_group_number {
        1..=2 => "P",
        3..=15 if monoclinic_is_c(space_group_number) => "C",
        3..=15 => "P",
        16..=74 => orthorhombic_centring(space_group_number),
        75..=142 if tetragonal_is_i(space_group_number) => "I",
        75..=142 => "P",
        143..=167 if is_rhombohedral(space_group_number) => "R",
        143..=194 => "P",
        _ if cubic_is_f(space_group_number) => "F",
        _ if cubic_is_i(space_group_number) => "I",
        _ => "P",
    })
}

fn monoclinic_is_c(number: u32) -> bool {
    matches!(number, 5 | 8 | 9 | 12 | 15)
}

fn is_rhombohedral(number: u32) -> bool {
    matches!(number, 146 | 148 | 155 | 160 | 161 | 166 | 167)
}

fn orthorhombic_centring(number: u32) -> &'static str {
    if matches!(number, 38..=41) {
        "A"
    } else if matches!(number, 20 | 21 | 35 | 36 | 37 | 63..=68) {
        "C"
    } else if matches!(number, 22 | 42 | 43 | 69 | 70) {
        "F"
    } else if matches!(number, 23 | 24 | 44 | 45 | 46 | 71..=74) {
        "I"
    } else {
        "P"
    }
}

fn tetragonal_is_i(number: u32) -> bool {
    matches!(
        number,
        79 | 80
            | 82
            | 87
            | 88
            | 97
            | 98
            | 107
            | 108
            | 109
            | 110
            | 119
            | 120
            | 121
            | 122
            | 139
            | 140
            | 141
            | 142
    )
}

fn cubic_is_f(number: u32) -> bool {
    matches!(
        number,
        196 | 202 | 203 | 209 | 210 | 216 | 219 | 225 | 226 | 227 | 228
    )
}

fn cubic_is_i(number: u32) -> bool {
    matches!(
        number,
        197 | 199 | 204 | 206 | 211 | 214 | 217 | 220 | 229 | 230
    )
}
