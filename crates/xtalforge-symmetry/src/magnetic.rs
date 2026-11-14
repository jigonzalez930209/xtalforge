//! Collinear magnetic space groups and moment symmetrisation (Tasks 2.8.2 and
//! 19.1.2).
//!
//! # Deviation from spglib golden data (Task 19.1.2)
//!
//! The spglib magnetic test suite is not redistributed with this repository,
//! so the tests use deterministic hand-built cases (ferromagnetic bcc Fe, a
//! checkerboard antiferromagnet on a simple-cubic supercell, a two-sublattice
//! ferrimagnet) and check orbit/stabiliser structure, operation counts and
//! moment symmetrisation instead of comparing UNI numbers with spglib golden
//! files. The backend is moyo; its non-magnetic datasets are audited against
//! spglib in `crates/xtalforge-symmetry/tests/moyo_audit.rs`.

use moyo::MoyoMagneticDataset;
use moyo::base::{AngleTolerance, Collinear, MagneticCell, RotationMagneticMomentAction};
use serde::Serialize;
use serde_json::Value;
use xtalforge_core::{Error, MagneticMoment, P1Cell, Structure};

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

/// Options for a magnetic symmetry search (Task 19.1.2).
#[derive(Debug, Clone, PartialEq)]
pub struct MagneticSymmetryOptions {
    /// Positional tolerance in Angstrom.
    pub symprec: f64,
    /// Moment tolerance in Bohr magnetons; `None` reuses `symprec`.
    ///
    /// The default is `Some(0.1)`, which absorbs typical DFT moment noise
    /// before the group is averaged.
    pub mag_symprec: Option<f64>,
    /// Rotate the standardised basis into the canonical orientation.
    pub rotate_basis: bool,
}

impl Default for MagneticSymmetryOptions {
    fn default() -> Self {
        Self {
            symprec: 1e-3,
            mag_symprec: Some(0.1),
            rotate_basis: false,
        }
    }
}

/// Result of symmetrising a collinear moment array (Task 19.1.2).
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct MagneticSymmetrization {
    /// UNI number of the detected magnetic space group.
    pub uni_number: u32,
    /// Number of magnetic operations in the input cell.
    pub operation_count: usize,
    /// Number of those operations that apply time reversal.
    pub time_reversal_count: usize,
    /// Symmetrised moments in P1 atom order.
    pub moments: Vec<f64>,
    /// Representative atom index for each input atom.
    pub orbits: Vec<usize>,
    /// Number of operations that leave each atom in place (site stabiliser
    /// order), in P1 atom order.
    pub stabilizers: Vec<usize>,
}

/// Detects the collinear magnetic space group of a structure.
///
/// Magnetic moments come from the typed site moment when present, otherwise
/// from the site property `magmom` (a number); a site without either is
/// treated as non-magnetic (`0`). A non-numeric property is also treated as
/// non-magnetic, matching the historical behaviour. The action on the moments
/// is axial.
///
/// # Errors
///
/// Returns [`Error::Symmetry`] when the search fails or the tolerance is
/// invalid.
pub fn detect_magnetic_symmetry(
    structure: &Structure,
    symprec: f64,
) -> Result<MagneticDataset, Error> {
    let p1 = structure.p1()?;
    let moments = lenient_moments(structure, &p1);
    let options = MagneticSymmetryOptions {
        symprec,
        mag_symprec: None,
        rotate_basis: false,
    };
    run_search(&p1, &moments, &options)
}

/// Detects the collinear magnetic space group of a structure with explicit
/// options (Task 19.1.2).
///
/// Unlike [`detect_magnetic_symmetry`], a non-collinear typed moment or a
/// non-numeric `magmom` is rejected: use this entry point when the input must
/// be explicitly collinear.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-collinear `magmom` and
/// [`Error::Symmetry`] when the search fails or an option is invalid.
pub fn detect_magnetic_symmetry_with_options(
    structure: &Structure,
    options: &MagneticSymmetryOptions,
) -> Result<MagneticDataset, Error> {
    let p1 = structure.p1()?;
    let moments = collinear_moments(structure, &p1)?;
    run_search(&p1, &moments, options)
}

/// Symmetrises a collinear moment array under the detected magnetic group.
///
/// `moments` are in the order of `structure.p1()?.atoms`. The array is
/// projected onto the invariant subspace of the magnetic operations: every
/// operation maps atom `i` to atom `j` and transforms the axial moment by the
/// determinant of its rotation and by time reversal, then the contributions
/// to each atom are averaged. Atoms belonging to one magnetic orbit therefore
/// carry the same symmetrised magnitude.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] when the array length differs from the P1
/// atom count, a moment is non-finite, or an option is invalid, and
/// [`Error::Symmetry`] when the search or the operation-to-atom mapping fails.
pub fn symmetrize_moments(
    structure: &Structure,
    moments: &[f64],
    options: &MagneticSymmetryOptions,
) -> Result<MagneticSymmetrization, Error> {
    validate_options(options)?;
    let p1 = structure.p1()?;
    if moments.len() != p1.atoms.len() {
        return Err(Error::InvalidInput(format!(
            "expected {} moments for the P1 cell, got {}",
            p1.atoms.len(),
            moments.len()
        )));
    }
    if !moments.iter().all(|value| value.is_finite()) {
        return Err(Error::InvalidInput(
            "magnetic moments must be finite".to_owned(),
        ));
    }
    let dataset = run_search(&p1, moments, options)?;
    let atoms = p1.atoms.len();
    let mut sums = vec![0.0_f64; atoms];
    let mut counts = vec![0_usize; atoms];
    let mut stabilizers = vec![0_usize; atoms];
    let tolerance = options.symprec.max(1e-6);
    for operation in &dataset.operations {
        let determinant = determinant(&operation.rotation)?;
        let time_reversal = if operation.time_reversal { -1 } else { 1 };
        let sign = f64::from(determinant * time_reversal);
        for (index, (atom, moment)) in p1.atoms.iter().zip(moments).enumerate() {
            let image = apply_operation(operation, atom.frac);
            let target = find_image(&p1, atom, image, tolerance)?;
            let Some(sum) = sums.get_mut(target) else {
                return Err(Error::Symmetry(format!(
                    "magnetic operation targets missing atom index {target}"
                )));
            };
            let Some(count) = counts.get_mut(target) else {
                return Err(Error::Symmetry(format!(
                    "magnetic operation targets missing atom index {target}"
                )));
            };
            *sum += sign * *moment;
            *count += 1;
            if target == index
                && let Some(stabilizer) = stabilizers.get_mut(index)
            {
                *stabilizer += 1;
            }
        }
    }
    let mut averaged = Vec::with_capacity(atoms);
    for (index, (sum, count)) in sums.iter().zip(&counts).enumerate() {
        if *count == 0 {
            return Err(Error::Symmetry(format!(
                "magnetic operations do not map atom {index} onto the cell"
            )));
        }
        averaged.push(sum / *count as f64);
    }
    Ok(MagneticSymmetrization {
        uni_number: dataset.uni_number,
        operation_count: dataset.operations.len(),
        time_reversal_count: dataset
            .operations
            .iter()
            .filter(|operation| operation.time_reversal)
            .count(),
        moments: averaged,
        orbits: dataset.orbits,
        stabilizers,
    })
}

/// Symmetrises the `magmom` site properties of a structure.
///
/// This is [`symmetrize_moments`] with the moments read from the site
/// properties. A missing property is `0`; a non-numeric property is rejected
/// because only collinear moments are supported here.
///
/// # Errors
///
/// Returns the errors of [`symmetrize_moments`] and [`Error::InvalidInput`]
/// for a non-collinear `magmom`.
pub fn symmetrize_structure_moments(
    structure: &Structure,
    options: &MagneticSymmetryOptions,
) -> Result<MagneticSymmetrization, Error> {
    let moments = magmom_values(structure)?;
    symmetrize_moments(structure, &moments, options)
}

/// Reads the collinear moments in P1 atom order.
///
/// The typed site moment wins when present; otherwise the `magmom` property is
/// used. A missing value is `0`. A non-collinear typed moment or a non-numeric
/// property is rejected.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-numeric `magmom` or a non-finite
/// number.
pub fn magmom_values(structure: &Structure) -> Result<Vec<f64>, Error> {
    let p1 = structure.p1()?;
    collinear_moments(structure, &p1)
}

/// Formats collinear moments as the space-separated value of a VASP `MAGMOM`
/// tag (four decimals, matching the Phase 18 writer).
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty array or a non-finite value.
pub fn format_magmom(moments: &[f64]) -> Result<String, Error> {
    if moments.is_empty() {
        return Err(Error::InvalidInput(
            "MAGMOM needs at least one moment".to_owned(),
        ));
    }
    if !moments.iter().all(|value| value.is_finite()) {
        return Err(Error::InvalidInput(
            "MAGMOM moments must be finite".to_owned(),
        ));
    }
    Ok(moments
        .iter()
        .map(|value| format!("{value:.4}"))
        .collect::<Vec<_>>()
        .join(" "))
}

fn validate_options(options: &MagneticSymmetryOptions) -> Result<(), Error> {
    if !options.symprec.is_finite() || options.symprec <= 0.0 {
        return Err(Error::Symmetry(format!(
            "symprec must be positive and finite, got {}",
            options.symprec
        )));
    }
    if let Some(mag_symprec) = options.mag_symprec
        && (!mag_symprec.is_finite() || mag_symprec <= 0.0)
    {
        return Err(Error::Symmetry(format!(
            "mag_symprec must be positive and finite, got {mag_symprec}"
        )));
    }
    Ok(())
}

fn run_search(
    p1: &P1Cell,
    moments: &[f64],
    options: &MagneticSymmetryOptions,
) -> Result<MagneticDataset, Error> {
    validate_options(options)?;
    if moments.len() != p1.atoms.len() {
        return Err(Error::InvalidInput(format!(
            "expected {} magnetic moments, got {}",
            p1.atoms.len(),
            moments.len()
        )));
    }
    let moyo_cell = to_moyo(p1)?;
    let magnetic_moments = moments.iter().copied().map(Collinear).collect();
    let magnetic_cell = MagneticCell::from_cell(moyo_cell.cell, magnetic_moments);
    let dataset = MoyoMagneticDataset::new(
        &magnetic_cell,
        options.symprec,
        AngleTolerance::Default,
        options.mag_symprec,
        RotationMagneticMomentAction::Axial,
        options.rotate_basis,
    )
    .map_err(|error| {
        Error::Symmetry(format!(
            "magnetic search failed at symprec={} A: {error}",
            options.symprec
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

fn lenient_moments(structure: &Structure, p1: &P1Cell) -> Vec<f64> {
    p1.atoms
        .iter()
        .map(|atom| {
            let site = structure.sites().get(atom.site_idx);
            if let Some(MagneticMoment::Collinear { magnitude, .. }) =
                site.and_then(|site| site.magnetic_moment())
            {
                return if magnitude.is_finite() {
                    magnitude
                } else {
                    0.0
                };
            }
            site.and_then(|site| site.props().get("magmom"))
                .and_then(Value::as_f64)
                .unwrap_or(0.0)
        })
        .collect()
}

fn collinear_moments(structure: &Structure, p1: &P1Cell) -> Result<Vec<f64>, Error> {
    p1.atoms
        .iter()
        .map(|atom| {
            let site = structure.sites().get(atom.site_idx);
            if let Some(moment) = site.and_then(|site| site.magnetic_moment()) {
                return match moment {
                    MagneticMoment::Collinear { magnitude, .. } => {
                        if magnitude.is_finite() {
                            Ok(magnitude)
                        } else {
                            Err(Error::InvalidInput(
                                "magnetic moments must be finite".to_owned(),
                            ))
                        }
                    }
                    MagneticMoment::NonCollinear { .. } => Err(Error::InvalidInput(format!(
                        "atom {} carries a non-collinear magnetic moment; only scalar moments are supported",
                        atom.label
                    ))),
                };
            }
            match site.and_then(|site| site.props().get("magmom")) {
                None | Some(Value::Null) => Ok(0.0),
                Some(Value::Number(number)) => number.as_f64().map_or_else(
                    || {
                        Err(Error::InvalidInput(
                            "a magnetic moment is not representable as f64".to_owned(),
                        ))
                    },
                    |value| {
                        if value.is_finite() {
                            Ok(value)
                        } else {
                            Err(Error::InvalidInput(
                                "magnetic moments must be finite".to_owned(),
                            ))
                        }
                    },
                ),
                Some(_) => Err(Error::InvalidInput(format!(
                    "atom {} carries a non-collinear magmom; only scalar moments are supported",
                    atom.label
                ))),
            }
        })
        .collect()
}

fn apply_operation(operation: &MagneticOperation, frac: [f64; 3]) -> [f64; 3] {
    let [[a, b, c], [d, e, f], [g, h, i]] = operation.rotation;
    [
        f64::from(a) * frac[0]
            + f64::from(b) * frac[1]
            + f64::from(c) * frac[2]
            + operation.translation[0],
        f64::from(d) * frac[0]
            + f64::from(e) * frac[1]
            + f64::from(f) * frac[2]
            + operation.translation[1],
        f64::from(g) * frac[0]
            + f64::from(h) * frac[1]
            + f64::from(i) * frac[2]
            + operation.translation[2],
    ]
}

fn find_image(
    p1: &P1Cell,
    atom: &xtalforge_core::P1Atom,
    image: [f64; 3],
    tolerance: f64,
) -> Result<usize, Error> {
    let mut best: Option<(usize, f64)> = None;
    for (index, candidate) in p1.atoms.iter().enumerate() {
        if candidate.element != atom.element || (candidate.occupancy - atom.occupancy).abs() > 1e-6
        {
            continue;
        }
        let distance = p1.lattice.distance(candidate.frac, image);
        if distance > tolerance {
            continue;
        }
        if best.is_none_or(|(_, current)| distance < current) {
            best = Some((index, distance));
        }
    }
    best.map(|(index, _)| index).ok_or_else(|| {
        Error::Symmetry(format!(
            "no image of atom {} at [{:.6}, {:.6}, {:.6}] within {tolerance} A",
            atom.label, image[0], image[1], image[2]
        ))
    })
}

fn determinant(rotation: &[[i32; 3]; 3]) -> Result<i32, Error> {
    let [[a, b, c], [d, e, f], [g, h, i]] = *rotation;
    let value = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    if value == 1 || value == -1 {
        Ok(value)
    } else {
        Err(Error::Symmetry(format!(
            "magnetic operation has a non-unimodular rotation (det={value})"
        )))
    }
}
