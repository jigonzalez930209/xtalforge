//! Tolerance ladders and symmetry diagnostics (Task 2.9).

use serde::Serialize;
use xtalforge_core::{Error, Structure};

use crate::dataset::{SymmetryOptions, detect_symmetry};
use crate::spacegroup::SpaceGroup;

/// A range of tolerances that gives the same detected setting.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ToleranceRange {
    /// Smallest tolerance in the range.
    pub min: f64,
    /// Largest tolerance in the range.
    pub max: f64,
    /// Space-group number.
    pub number: u32,
    /// Hall number.
    pub hall_number: u32,
    /// H-M symbol.
    pub hm_symbol: String,
}

/// Detection result at one tolerance.
#[derive(Debug, Clone, PartialEq)]
pub struct ToleranceSample {
    /// Tolerance in Angstrom.
    pub symprec: f64,
    /// Space-group number.
    pub number: u32,
    /// Hall number.
    pub hall_number: u32,
    /// H-M symbol.
    pub hm_symbol: String,
}

/// Default tolerance ladder from 1e-5 to 0.5 Angstrom.
pub const DEFAULT_LADDER: [f64; 10] = [1e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 0.1, 0.2, 0.5];

/// Detects the structure at every tolerance of the ladder.
///
/// # Errors
///
/// Returns the first detection error.
pub fn tolerance_samples(
    structure: &Structure,
    tolerances: &[f64],
) -> Result<Vec<ToleranceSample>, Error> {
    let mut samples = Vec::with_capacity(tolerances.len());
    for symprec in tolerances {
        let dataset = detect_symmetry(
            structure,
            SymmetryOptions {
                symprec: *symprec,
                ..SymmetryOptions::default()
            },
        )?;
        samples.push(ToleranceSample {
            symprec: *symprec,
            number: dataset.number,
            hall_number: dataset.hall_number,
            hm_symbol: dataset.hm_symbol,
        });
    }
    Ok(samples)
}

/// Groups consecutive samples that share a Hall number into ranges.
///
/// # Errors
///
/// Returns the errors of [`tolerance_samples`].
pub fn tolerance_ladder(
    structure: &Structure,
    tolerances: &[f64],
) -> Result<Vec<ToleranceRange>, Error> {
    let samples = tolerance_samples(structure, tolerances)?;
    let mut ranges: Vec<ToleranceRange> = Vec::new();
    for sample in samples {
        if let Some(last) = ranges.last_mut()
            && last.hall_number == sample.hall_number
        {
            last.max = sample.symprec;
            continue;
        }
        ranges.push(ToleranceRange {
            min: sample.symprec,
            max: sample.symprec,
            number: sample.number,
            hall_number: sample.hall_number,
            hm_symbol: sample.hm_symbol,
        });
    }
    Ok(ranges)
}

/// Deviation of one atom from its orbit average.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct AtomDeviation {
    /// Atom index in the input cell.
    pub atom: usize,
    /// Element symbol.
    pub element: String,
    /// Site label.
    pub label: String,
    /// Cartesian displacement in Angstrom.
    pub displacement: f64,
}

/// Per-atom deviations from a requested space group.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SymmetryDiagnostics {
    /// Largest displacement.
    pub max_displacement: f64,
    /// Mean displacement.
    pub mean_displacement: f64,
    /// Per-atom records, sorted by decreasing displacement.
    pub atoms: Vec<AtomDeviation>,
}

/// Measures how far each atom is from the orbit average of a group.
///
/// # Errors
///
/// Returns [`Error::Symmetry`] for an empty structure or invalid group.
pub fn diagnose(
    structure: &Structure,
    space_group: &SpaceGroup,
) -> Result<SymmetryDiagnostics, Error> {
    let p1 = structure.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::Symmetry(
            "symmetry diagnostics need at least one atom".to_owned(),
        ));
    }
    if space_group.operations.is_empty() {
        return Err(Error::Symmetry(
            "symmetry diagnostics need a group with operations".to_owned(),
        ));
    }
    let mut atoms = Vec::with_capacity(p1.atoms.len());
    let mut sum = 0.0;
    let mut max = 0.0_f64;
    for (index, atom) in p1.atoms.iter().enumerate() {
        let mut deviation = f64::INFINITY;
        for operation in &space_group.operations {
            let image = operation.apply(atom.frac);
            for other in &p1.atoms {
                if other.element != atom.element {
                    continue;
                }
                deviation = deviation.min(p1.lattice.distance(image, other.frac));
            }
        }
        let deviation = if deviation.is_finite() {
            deviation
        } else {
            0.0
        };
        sum += deviation;
        max = max.max(deviation);
        atoms.push(AtomDeviation {
            atom: index,
            element: atom.element.clone(),
            label: atom.label.clone(),
            displacement: deviation,
        });
    }

    Ok(SymmetryDiagnostics {
        max_displacement: max,
        mean_displacement: sum / p1.atoms.len() as f64,
        atoms,
    })
}
