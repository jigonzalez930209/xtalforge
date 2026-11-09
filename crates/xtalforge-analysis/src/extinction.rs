//! Space-group screening from systematic absences (Task 7.10.1).
//!
//! The automatic run indexes the metric families first. The indexed
//! reflections then decide which space groups are compatible with what the
//! pattern actually shows: a reflection that a candidate group makes
//! systematically absent cannot produce an observed line.
//!
//! # Method
//!
//! For one candidate group the structure factor of an atom on a general
//! position is
//!
//! ```text
//! F(hkl) = sum over operations of exp(2 pi i (h R r + h t))
//! ```
//!
//! which is zero for every `r` exactly when the reflection is systematically
//! absent. Two generic trial points evaluate that sum; a candidate is
//! rejected for the reflections whose `|F|` is numerically zero at both.
//!
//! # Approximations
//!
//! * The screen is a reflection-condition test, not the reference's Le Bail
//!   class fit. It cannot rank two groups that share the same absences by
//!   fit quality; they are ordered by the number of operations, so the
//!   higher-symmetry group comes first.
//! * Only the crystal system of the indexed metric family is searched (plus
//!   trigonal groups for a hexagonal cell). A lower-symmetry subgroup in a
//!   smaller cell is not considered.
//! * Observed peaks that the indexing could not explain at all are not in
//!   the list, so they cannot make a candidate look extinct. A peak assigned
//!   to the wrong reflection can.
//! * When no candidate explains every indexed reflection, the groups with
//!   the fewest unexplained reflections are returned and the count is
//!   reported; an empty list means no space group could be screened.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, SymOp};
use xtalforge_symmetry::spacegroup::{SpaceGroup, SpaceGroupFilter, catalogue};

use crate::indexing::LatticeSystem;

/// Maximum number of candidates reported per cell.
pub const MAX_CANDIDATES: usize = 8;

/// `|F|` below this is a systematic absence at the trial points.
const ABSENT: f64 = 1e-6;

/// Two generic fractional points used to evaluate the structure factor.
const TRIAL_POINTS: [[f64; 3]; 2] = [[0.137, 0.281, 0.419], [0.217, 0.353, 0.071]];

/// One space group compatible with the indexed reflections.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SpaceGroupMatch {
    /// ITA number (1 to 230).
    pub number: u32,
    /// Hermann-Mauguin short symbol.
    pub symbol: String,
    /// Centering letter (`P`, `A`, `B`, `C`, `I`, `R` or `F`).
    pub centering: char,
    /// Geometric crystal class, for example `4/mmm`.
    pub point_group: String,
    /// Number of operations in the conventional cell.
    pub order: usize,
    /// Indexed reflections the group makes systematically absent.
    pub unexplained_peaks: usize,
}

/// Whether a group makes the reflection systematically absent.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty operation list, which cannot
/// describe a space group.
pub fn is_systematically_absent(operations: &[SymOp], hkl: [i32; 3]) -> Result<bool, Error> {
    if operations.is_empty() {
        return Err(Error::InvalidInput(
            "a space group needs at least one operation".to_owned(),
        ));
    }
    for point in TRIAL_POINTS {
        let mut real = 0.0;
        let mut imaginary = 0.0;
        for operation in operations {
            let image = operation.apply(point);
            let phase = std::f64::consts::TAU
                * (f64::from(hkl[0]) * image[0]
                    + f64::from(hkl[1]) * image[1]
                    + f64::from(hkl[2]) * image[2]);
            real += phase.cos();
            imaginary += phase.sin();
        }
        let scale = operations.len() as f64;
        let magnitude = (real * real + imaginary * imaginary).sqrt() / scale;
        if magnitude > ABSENT {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Crystal-system names searched for an indexed metric family.
fn systems_of(system: LatticeSystem) -> &'static [&'static str] {
    match system {
        LatticeSystem::Cubic => &["cubic"],
        LatticeSystem::Tetragonal => &["tetragonal"],
        LatticeSystem::Hexagonal => &["hexagonal", "trigonal"],
        LatticeSystem::Orthorhombic => &["orthorhombic"],
        LatticeSystem::Monoclinic => &["monoclinic"],
    }
}

/// Screens the space groups of a metric family against the indexed hkl.
///
/// # Errors
///
/// Returns [`Error::Symmetry`] if the space-group catalogue cannot be read,
/// and the errors of [`is_systematically_absent`].
pub fn space_group_candidates(
    system: LatticeSystem,
    hkl: &[[i32; 3]],
) -> Result<Vec<SpaceGroupMatch>, Error> {
    let mut groups: BTreeMap<u32, SpaceGroup> = BTreeMap::new();
    for name in systems_of(system) {
        for group in catalogue(SpaceGroupFilter {
            crystal_system: Some(name),
            ..SpaceGroupFilter::default()
        })? {
            groups.entry(group.number).or_insert(group);
        }
    }
    let mut matches = Vec::new();
    for group in groups.into_values() {
        let mut unexplained = 0usize;
        for reflection in hkl {
            if *reflection == [0, 0, 0] {
                continue;
            }
            if is_systematically_absent(&group.operations, *reflection)? {
                unexplained += 1;
            }
        }
        matches.push(SpaceGroupMatch {
            number: group.number,
            symbol: group.hm_short,
            centering: group.centering,
            point_group: group.point_group,
            order: group.operations.len(),
            unexplained_peaks: unexplained,
        });
    }
    let fewest = matches
        .iter()
        .map(|entry| entry.unexplained_peaks)
        .min()
        .unwrap_or(0);
    matches.retain(|entry| entry.unexplained_peaks == fewest);
    matches.sort_by(|first, second| {
        second
            .order
            .cmp(&first.order)
            .then_with(|| first.number.cmp(&second.number))
    });
    matches.truncate(MAX_CANDIDATES);
    Ok(matches)
}
