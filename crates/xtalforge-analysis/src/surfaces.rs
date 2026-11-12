#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Symmetry-unique Miller indices and slab terminations (Tasks 9.4.1-9.4.2).

use std::collections::BTreeSet;

use serde::{Deserialize, Serialize};

use xtalforge_core::{Error, Lattice, Site, Slab, Structure, StructureBuilder, make_slab};
use xtalforge_symmetry::{SymmetryOptions, detect_symmetry};

use crate::bvs::{BvsOptions, oxidation_states};

/// One symmetry-unique Miller index.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MillerIndex {
    /// Representative `(h k l)`, reduced by the greatest common divisor.
    pub hkl: [i32; 3],
    /// Interplanar spacing in Angstrom.
    pub d_spacing: f64,
    /// Number of equivalent directions.
    pub multiplicity: usize,
    /// All equivalent directions, including the representative.
    pub family: Vec<[i32; 3]>,
}

/// Enumerates the symmetry-unique `(h k l)` up to `max_index`.
///
/// The point group is taken from a symmetry search of the structure and
/// inversion is always included (Laue class), as in pymatgen's
/// `get_symmetrically_distinct_miller_indices`.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a negative maximum index, or the errors
/// of [`detect_symmetry`] and [`Lattice::d_spacing`].
pub fn unique_miller_indices(
    structure: &Structure,
    max_index: i32,
    symprec: f64,
) -> Result<Vec<MillerIndex>, Error> {
    if max_index < 1 {
        return Err(Error::InvalidInput(
            "the Miller index bound must be positive".to_owned(),
        ));
    }
    let dataset = detect_symmetry(
        structure,
        SymmetryOptions {
            symprec,
            ..SymmetryOptions::default()
        },
    )?;
    let rotations: Vec<[[i32; 3]; 3]> = dataset
        .operations
        .iter()
        .map(|operation| operation.rotation())
        .collect();
    let mut seen: BTreeSet<[i32; 3]> = BTreeSet::new();
    let mut indices: Vec<MillerIndex> = Vec::new();
    for h in -max_index..=max_index {
        for k in -max_index..=max_index {
            for l in -max_index..=max_index {
                if h == 0 && k == 0 && l == 0 {
                    continue;
                }
                let gcd = gcd(gcd(h.abs(), k.abs()), l.abs()).max(1);
                if gcd > 1 {
                    continue;
                }
                let mut family = BTreeSet::new();
                for rotation in &rotations {
                    for sign in [1, -1] {
                        let transformed =
                            multiply_rotation(*rotation, [h * sign, k * sign, l * sign]);
                        family.insert(transformed);
                    }
                }
                let Some(canonical) = family
                    .iter()
                    .copied()
                    .filter(|vector| canonical_direction(*vector) == *vector)
                    .min_by_key(|vector| {
                        let positive = vector.iter().filter(|value| **value > 0).count() as i32;
                        (
                            -positive,
                            -vector[0].abs(),
                            -vector[1].abs(),
                            -vector[2].abs(),
                            *vector,
                        )
                    })
                else {
                    continue;
                };
                if seen.contains(&canonical) {
                    continue;
                }
                seen.insert(canonical);
                let d_spacing = structure.lattice().d_spacing(canonical.map(f64::from))?;
                indices.push(MillerIndex {
                    hkl: canonical,
                    d_spacing,
                    multiplicity: family.len(),
                    family: family.into_iter().collect(),
                });
            }
        }
    }
    indices.sort_by(|first, second| {
        second
            .d_spacing
            .partial_cmp(&first.d_spacing)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| first.hkl.cmp(&second.hkl))
    });
    Ok(indices)
}

/// One distinct termination of a slab.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Termination {
    /// Index in the report.
    pub index: usize,
    /// Shift passed to [`make_slab`].
    pub shift: f64,
    /// Number of atoms in the slab.
    pub n_atoms: usize,
    /// Slab thickness in Angstrom.
    pub thickness: f64,
    /// Net dipole moment along the surface normal in e Angstrom.
    pub dipole: f64,
    /// Whether the dipole is non-zero (a polar surface).
    pub polar: bool,
    /// The generated (possibly symmetric) slab.
    pub structure: Structure,
}

/// All terminations of one `(h k l)` surface.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct TerminationReport {
    /// Reduced Miller index.
    pub hkl: [i32; 3],
    /// Number of distinct cuts found.
    pub n_terminations: usize,
    /// Distinct terminations, sorted by shift.
    pub terminations: Vec<Termination>,
}

/// Generates every distinct termination of a `(h k l)` slab.
///
/// The candidate cuts are the midpoints between consecutive atomic layers in
/// the surface basis; duplicates after the cut are removed.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a bad layer count or vacuum, or the
/// errors of [`make_slab`] and [`oxidation_states`].
pub fn slab_terminations(
    structure: &Structure,
    hkl: [i32; 3],
    layers: u32,
    vacuum: f64,
    symmetric: bool,
) -> Result<TerminationReport, Error> {
    if layers == 0 {
        return Err(Error::InvalidInput(
            "a termination study needs at least one layer".to_owned(),
        ));
    }
    let basis = xtalforge_core::surface_basis(hkl)?;
    let basis_float = basis.map(|row| row.map(f64::from));
    let inverse = Lattice::new(basis_float)?.inverse();
    let p1 = structure.p1()?;
    let mut heights: Vec<f64> = p1
        .atoms
        .iter()
        .map(|atom| {
            let shifted = [
                atom.frac[0] * inverse[0][0]
                    + atom.frac[1] * inverse[1][0]
                    + atom.frac[2] * inverse[2][0],
                atom.frac[0] * inverse[0][1]
                    + atom.frac[1] * inverse[1][1]
                    + atom.frac[2] * inverse[2][1],
                atom.frac[0] * inverse[0][2]
                    + atom.frac[1] * inverse[1][2]
                    + atom.frac[2] * inverse[2][2],
            ];
            shifted[2].rem_euclid(1.0)
        })
        .collect();
    heights.sort_by(|first, second| {
        first
            .partial_cmp(second)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    heights.dedup_by(|first, second| (*first - *second).abs() < 1e-4);
    if heights.is_empty() {
        return Err(Error::InvalidInput(
            "the structure has no layers in this direction".to_owned(),
        ));
    }
    let mut shifts = Vec::with_capacity(heights.len());
    for (index, height) in heights.iter().enumerate() {
        let next = heights
            .get(index + 1)
            .copied()
            .unwrap_or_else(|| heights.first().copied().unwrap_or(0.0) + 1.0);
        shifts.push(((height + next) / 2.0).rem_euclid(1.0));
    }
    let mut terminations: Vec<Termination> = Vec::new();
    let mut signatures: BTreeSet<String> = BTreeSet::new();
    for shift in shifts {
        let slab = make_slab(structure, hkl, layers, vacuum, shift)?;
        let slab = if symmetric { mirror_slab(&slab)? } else { slab };
        let signature = slab_signature(&slab.structure);
        if !signatures.insert(signature) {
            continue;
        }
        let dipole = slab_dipole(&slab)?;
        let index = terminations.len();
        terminations.push(Termination {
            index,
            shift,
            n_atoms: slab.structure.sites().len(),
            thickness: slab.thickness,
            dipole,
            polar: dipole.abs() > 1e-4,
            structure: slab.structure,
        });
    }
    Ok(TerminationReport {
        hkl: hkl.map(|value| {
            let divisor = gcd(gcd(hkl[0].abs(), hkl[1].abs()), hkl[2].abs()).max(1);
            value / divisor
        }),
        n_terminations: terminations.len(),
        terminations,
    })
}

/// Builds the thinnest slab that reaches `thickness` Angstrom.
///
/// # Errors
///
/// Returns the errors of [`slab_terminations`] with one candidate cut.
pub fn minimal_slab(
    structure: &Structure,
    hkl: [i32; 3],
    thickness: f64,
    vacuum: f64,
) -> Result<Slab, Error> {
    if !thickness.is_finite() || thickness <= 0.0 {
        return Err(Error::InvalidInput(
            "the slab thickness must be positive and finite".to_owned(),
        ));
    }
    let spacing = structure.lattice().d_spacing(hkl.map(f64::from))?;
    let layers = (thickness / spacing).ceil().max(1.0) as u32;
    make_slab(structure, hkl, layers, vacuum, 0.0)
}

fn mirror_slab(slab: &Slab) -> Result<Slab, Error> {
    let p1 = slab.structure.p1()?;
    let min = p1
        .atoms
        .iter()
        .map(|atom| atom.frac[2])
        .fold(f64::INFINITY, f64::min);
    let max = p1
        .atoms
        .iter()
        .map(|atom| atom.frac[2])
        .fold(f64::NEG_INFINITY, f64::max);
    let middle = (min + max) / 2.0;
    let mut builder = StructureBuilder::new(slab.structure.name(), p1.lattice);
    for atom in &p1.atoms {
        builder = builder.site(
            Site::new(atom.label.clone(), &atom.element, atom.frac)?
                .with_occupancy(atom.occupancy)?,
        );
    }
    for atom in &p1.atoms {
        let mirrored = 2.0 * middle - atom.frac[2];
        let frac = [atom.frac[0], atom.frac[1], mirrored];
        builder = builder.site(
            Site::new(
                format!("{}_m", atom.label),
                &atom.element,
                Lattice::wrap(frac),
            )?
            .with_occupancy(atom.occupancy)?,
        );
    }
    let structure = builder.build()?;
    let mut output = slab.clone();
    output.structure = structure;
    Ok(output)
}

fn slab_dipole(slab: &Slab) -> Result<f64, Error> {
    let oxidation = oxidation_states(&slab.structure, &BvsOptions::default());
    let charges = match oxidation {
        Ok(report) => report
            .sites
            .iter()
            .map(|site| f64::from(site.state))
            .collect::<Vec<_>>(),
        Err(_) => slab
            .structure
            .p1()?
            .atoms
            .iter()
            .map(|atom| {
                slab.structure
                    .sites()
                    .iter()
                    .find(|site| site.label() == atom.label)
                    .and_then(Site::charge)
                    .unwrap_or(0.0)
            })
            .collect(),
    };
    let p1 = slab.structure.p1()?;
    let normal = slab.structure.lattice().plane_normal([0.0, 0.0, 1.0])?;
    let mut dipole = 0.0;
    for (atom, charge) in p1.atoms.iter().zip(charges) {
        let cart = p1.lattice.to_cart(atom.frac);
        dipole += charge * (cart[0] * normal[0] + cart[1] * normal[1] + cart[2] * normal[2]);
    }
    Ok(dipole)
}

fn slab_signature(structure: &Structure) -> String {
    let mut atoms: Vec<(u64, String, u64, u64)> = structure
        .sites()
        .iter()
        .map(|site| {
            (
                (site.frac()[2] * 1000.0).round() as u64,
                site.element().to_owned(),
                (site.frac()[0] * 1000.0).round() as u64,
                (site.frac()[1] * 1000.0).round() as u64,
            )
        })
        .collect();
    atoms.sort();
    atoms
        .into_iter()
        .map(|(z, element, x, y)| format!("{element}:{z}:{x}:{y}"))
        .collect::<Vec<_>>()
        .join("|")
}

fn canonical_direction(vector: [i32; 3]) -> [i32; 3] {
    for value in vector {
        if value > 0 {
            return vector;
        }
        if value < 0 {
            return [-vector[0], -vector[1], -vector[2]];
        }
    }
    vector
}

fn multiply_rotation(rotation: [[i32; 3]; 3], hkl: [i32; 3]) -> [i32; 3] {
    [
        rotation[0][0] * hkl[0] + rotation[0][1] * hkl[1] + rotation[0][2] * hkl[2],
        rotation[1][0] * hkl[0] + rotation[1][1] * hkl[1] + rotation[1][2] * hkl[2],
        rotation[2][0] * hkl[0] + rotation[2][1] * hkl[1] + rotation[2][2] * hkl[2],
    ]
}

fn gcd(first: i32, second: i32) -> i32 {
    let mut a = first;
    let mut b = second;
    while b != 0 {
        let remainder = a % b;
        a = b;
        b = remainder;
    }
    a.abs()
}
