//! Collinear magnetic ordering enumeration (Task 19.1.4).
//!
//! The input carries one collinear moment per site, either as the typed site
//! moment or as the `magmom` property. Its absolute value is the sublattice
//! amplitude; the module ignores the input sign and enumerates every sign
//! decoration of the magnetic atoms in a supercell, fixing the first magnetic
//! atom positive to remove the global spin flip. Each
//! decoration is classified as ferro-, antiferro- or ferrimagnetic from its net
//! moment and ranked by the detected magnetic space group: more magnetic
//! operations first, then kind (ferro, antiferro, ferri), UNI number and moment
//! vector.
//!
//! The unit-cell and enumeration limits bound the work: at most
//! [`MagneticOrderingOptions::max_candidates`] sign patterns are evaluated,
//! chosen by increasing Hamming weight, and at most
//! [`MagneticOrderingOptions::max_orderings`] results are returned. A pattern
//! whose magnetic group cannot be determined is still returned with
//! `uni_number: None` and ranked last.

use std::cmp::Ordering;

use serde::{Deserialize, Serialize};
use serde_json::Value;

use xtalforge_core::{
    Error, MagneticAxis, MagneticMoment, P1Cell, Site, Structure, StructureBuilder,
};
use xtalforge_symmetry::{
    MagneticSymmetryOptions, detect_magnetic_symmetry_with_options, format_magmom,
};

/// Ferro-, antiferro- or ferrimagnetic character of a collinear ordering.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum MagneticOrderingKind {
    /// Every nonzero moment has the same sign.
    Ferromagnetic,
    /// Both signs appear and the cell moment cancels.
    Antiferromagnetic,
    /// Both signs appear and the cell moment does not cancel.
    Ferrimagnetic,
}

/// One collinear magnetic ordering of a supercell.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MagneticOrdering {
    /// Ordering class derived from the moment signs.
    pub kind: MagneticOrderingKind,
    /// P1 supercell with the `magmom` site property set and a collinear typed
    /// moment along `c`.
    pub structure: Structure,
    /// Moments in P1 atom order.
    pub moments: Vec<f64>,
    /// Sum of the moments in Bohr magnetons.
    pub net_moment: f64,
    /// UNI number of the detected magnetic group; `None` when the search
    /// failed for this decoration.
    pub uni_number: Option<u32>,
    /// Number of magnetic operations in the supercell.
    pub operation_count: usize,
    /// Value of a VASP `MAGMOM` tag, formatted like the Phase 18 writer.
    pub magmom: String,
}

/// Options for [`magnetic_orderings`] and
/// [`magnetic_orderings_in_supercell`].
#[derive(Debug, Clone, PartialEq)]
pub struct MagneticOrderingOptions {
    /// Supercell repeats used by [`magnetic_orderings`].
    ///
    /// [`magnetic_orderings_in_supercell`] uses the supplied cell as-is and
    /// ignores this field.
    pub supercell: [u32; 3],
    /// Maximum number of orderings returned.
    pub max_orderings: usize,
    /// Maximum number of sign patterns evaluated.
    pub max_candidates: usize,
    /// Magnetic symmetry search options used for ranking.
    pub symmetry: MagneticSymmetryOptions,
}

impl Default for MagneticOrderingOptions {
    fn default() -> Self {
        Self {
            supercell: [1, 1, 1],
            max_orderings: 64,
            max_candidates: 256,
            symmetry: MagneticSymmetryOptions::default(),
        }
    }
}

/// Enumerates and ranks collinear magnetic orderings of a supercell.
///
/// The moment amplitude of each magnetic atom is `|magmom|` (or the magnitude
/// of the typed collinear site moment); atoms without either are non-magnetic.
/// The input sign is not used, so a ferro- and an antiferromagnetic
/// description of the same cell produce the same candidate set.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a zero supercell, a zero limit, a
/// non-collinear moment or a structure without magnetic atoms, and propagates
/// the supercell and P1 expansion errors.
pub fn magnetic_orderings(
    structure: &Structure,
    options: &MagneticOrderingOptions,
) -> Result<Vec<MagneticOrdering>, Error> {
    if options.supercell.contains(&0) {
        return Err(Error::InvalidInput(
            "magnetic ordering supercell repeats must be positive".to_owned(),
        ));
    }
    let supercell = structure.supercell(
        options.supercell[0],
        options.supercell[1],
        options.supercell[2],
    )?;
    enumerate_in_cell(&supercell, options)
}

/// Enumerates and ranks collinear magnetic orderings of a supplied supercell.
///
/// A non-P1 input is expanded to P1 first. The `supercell` field of `options`
/// is ignored; use it to record which reads produced the cell.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a zero limit, a non-collinear `magmom`
/// or a structure without magnetic atoms, and propagates the P1 expansion
/// errors.
pub fn magnetic_orderings_in_supercell(
    supercell: &Structure,
    options: &MagneticOrderingOptions,
) -> Result<Vec<MagneticOrdering>, Error> {
    enumerate_in_cell(supercell, options)
}

fn enumerate_in_cell(
    supercell: &Structure,
    options: &MagneticOrderingOptions,
) -> Result<Vec<MagneticOrdering>, Error> {
    if options.max_orderings == 0 {
        return Err(Error::InvalidInput(
            "max_orderings must be positive".to_owned(),
        ));
    }
    if options.max_candidates == 0 {
        return Err(Error::InvalidInput(
            "max_candidates must be positive".to_owned(),
        ));
    }
    let p1 = supercell.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "a magnetic ordering needs at least one atom".to_owned(),
        ));
    }
    let amplitudes = magnetic_amplitudes(supercell, &p1)?;
    let magnetic: Vec<usize> = amplitudes
        .iter()
        .enumerate()
        .filter_map(|(index, amplitude)| (*amplitude > 1e-9).then_some(index))
        .collect();
    if magnetic.is_empty() {
        return Err(Error::InvalidInput(
            "the structure has no magnetic atoms (nonzero magmom)".to_owned(),
        ));
    }
    let patterns = sign_patterns(magnetic.len() - 1, options.max_candidates);
    let mut results = Vec::with_capacity(patterns.len());
    for reversed in &patterns {
        let moments = decorate(&amplitudes, &magnetic, reversed);
        let kind = classify(&moments);
        let net_moment = clean_sum(&moments);
        let candidate = with_moments(&p1, &moments)?;
        let group = detect_magnetic_symmetry_with_options(&candidate, &options.symmetry).ok();
        let (uni_number, operation_count) = group.map_or((None, 0), |dataset| {
            (Some(dataset.uni_number), dataset.operations.len())
        });
        let magmom = format_magmom(&moments)?;
        results.push(MagneticOrdering {
            kind,
            structure: candidate,
            moments,
            net_moment,
            uni_number,
            operation_count,
            magmom,
        });
    }
    results.sort_by(compare_orderings);
    results.truncate(options.max_orderings);
    Ok(results)
}

fn magnetic_amplitudes(structure: &Structure, p1: &P1Cell) -> Result<Vec<f64>, Error> {
    p1.atoms
        .iter()
        .map(|atom| {
            let site = structure.sites().get(atom.site_idx);
            if let Some(moment) = site.and_then(|site| site.magnetic_moment()) {
                return match moment {
                    MagneticMoment::Collinear { magnitude, .. } => {
                        if magnitude.is_finite() {
                            Ok(magnitude.abs())
                        } else {
                            Err(Error::InvalidInput(
                                "magnetic moments must be finite".to_owned(),
                            ))
                        }
                    }
                    MagneticMoment::NonCollinear { .. } => Err(Error::InvalidInput(format!(
                        "atom {} carries a non-collinear magnetic moment; only collinear ordering is supported",
                        atom.label
                    ))),
                };
            }
            match site.and_then(|site| site.props().get("magmom")) {
                None | Some(Value::Null) => Ok(0.0),
                Some(Value::Number(number)) => {
                    let value = number.as_f64().ok_or_else(|| {
                        Error::InvalidInput(
                            "a magnetic moment is not representable as f64".to_owned(),
                        )
                    })?;
                    if !value.is_finite() {
                        return Err(Error::InvalidInput(
                            "magnetic moments must be finite".to_owned(),
                        ));
                    }
                    Ok(value.abs())
                }
                Some(_) => Err(Error::InvalidInput(format!(
                    "atom {} carries a non-collinear magmom; only collinear ordering is supported",
                    atom.label
                ))),
            }
        })
        .collect()
}

/// Sign patterns over `reversals` atoms after fixing the first magnetic atom
/// positive, in increasing Hamming weight and lexicographic order.
fn sign_patterns(reversals: usize, limit: usize) -> Vec<Vec<usize>> {
    let mut output = Vec::new();
    for weight in 0..=reversals {
        if output.len() >= limit {
            break;
        }
        combinations(reversals, weight, 0, &mut Vec::new(), &mut output, limit);
    }
    output
}

fn combinations(
    n: usize,
    weight: usize,
    start: usize,
    current: &mut Vec<usize>,
    output: &mut Vec<Vec<usize>>,
    limit: usize,
) {
    if output.len() >= limit {
        return;
    }
    if current.len() == weight {
        output.push(current.clone());
        return;
    }
    let remaining = weight - current.len();
    if remaining > n.saturating_sub(start) {
        return;
    }
    for index in start..n {
        if output.len() >= limit {
            return;
        }
        if n - index < remaining {
            break;
        }
        current.push(index);
        combinations(n, weight, index + 1, current, output, limit);
        current.pop();
    }
}

fn decorate(amplitudes: &[f64], magnetic: &[usize], reversed: &[usize]) -> Vec<f64> {
    let mut moments = vec![0.0; amplitudes.len()];
    for (position, atom) in magnetic.iter().enumerate() {
        let Some(moment) = amplitudes.get(*atom) else {
            continue;
        };
        // Position 0 is the fixed-positive atom; `reversed` indexes the others.
        let sign = if position > 0 && reversed.contains(&(position - 1)) {
            -1.0
        } else {
            1.0
        };
        if let Some(slot) = moments.get_mut(*atom) {
            *slot = sign * moment;
        }
    }
    moments
}

fn classify(moments: &[f64]) -> MagneticOrderingKind {
    let mut positive = false;
    let mut negative = false;
    let mut scale = 0.0;
    for moment in moments {
        if *moment > 0.0 {
            positive = true;
        } else if *moment < 0.0 {
            negative = true;
        }
        scale += moment.abs();
    }
    let net = clean_sum(moments).abs();
    if !negative || !positive {
        MagneticOrderingKind::Ferromagnetic
    } else if net <= 1e-8 * scale.max(1.0) {
        MagneticOrderingKind::Antiferromagnetic
    } else {
        MagneticOrderingKind::Ferrimagnetic
    }
}

fn clean_sum(moments: &[f64]) -> f64 {
    let sum: f64 = moments.iter().sum();
    if sum.abs() <= 1e-12 { 0.0 } else { sum }
}

fn with_moments(p1: &P1Cell, moments: &[f64]) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new("magnetic-ordering", p1.lattice);
    for (index, (atom, moment)) in p1.atoms.iter().zip(moments).enumerate() {
        let label = format!("{}{}", atom.element, index + 1);
        let mut site = Site::new(label, &atom.element, atom.frac)?
            .with_occupancy(atom.occupancy)?
            .with_prop("magmom", serde_json::json!(moment));
        if moment.abs() > 1e-9 {
            site = site.with_collinear_moment(MagneticAxis::C, *moment)?;
        }
        builder = builder.site(site);
    }
    builder.build()
}

fn compare_orderings(first: &MagneticOrdering, second: &MagneticOrdering) -> Ordering {
    second
        .operation_count
        .cmp(&first.operation_count)
        .then_with(|| first.kind.cmp(&second.kind))
        .then_with(|| first.uni_number.cmp(&second.uni_number))
        .then_with(|| compare_moments(&first.moments, &second.moments))
}

fn compare_moments(first: &[f64], second: &[f64]) -> Ordering {
    for (left, right) in first.iter().zip(second) {
        let order = left.total_cmp(right);
        if order != Ordering::Equal {
            return order;
        }
    }
    first.len().cmp(&second.len())
}
