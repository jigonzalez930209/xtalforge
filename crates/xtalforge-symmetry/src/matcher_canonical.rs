//! Standardized primitive cells for the structure matcher (Task 17.3.3).
//!
//! Standardization removes the cell choice but not every orientation: two
//! inputs can reach primitive cells related by a point-group operation of the
//! same lattice. This module standardizes through the conventional setting and
//! then picks a deterministic representative among the lattice automorphisms,
//! so the search always compares decorations in the same frame.

use xtalforge_core::{Error, Lattice, P1Cell, Site, Structure, StructureBuilder};

use crate::dataset::SymmetryOptions;
use crate::matcher::{MatcherOptions, SpeciesComparator};
use crate::reduction::{Reduction, reduce_structure};
use crate::standardize::{primitive_cell, standardize};

/// Quantized decoration used to choose a canonical orientation.
type AtomSignature = Vec<(String, i64, i64, i64)>;

/// Chosen automorphism and the decoration its inverse produces.
type CanonicalChoice = (AtomSignature, [[f64; 3]; 3]);

/// Relative error allowed when a matrix is tested as a lattice automorphism.
const AUTOMORPHISM_TOLERANCE: f64 = 1e-6;

/// Reduces a structure to its standardized primitive cell.
///
/// Falls back to the Niggli-reduced input cell when standardization does not
/// produce a non-empty cell, which keeps the matcher usable on structures that
/// the symmetry search cannot handle.
pub(crate) fn prepare(structure: &Structure, options: &MatcherOptions) -> Result<P1Cell, Error> {
    if let Ok(cell) = standardized_primitive_cell(structure, options.comparator)
        && !cell.atoms.is_empty()
    {
        return Ok(cell);
    }
    let reduced = reduce_structure(structure, Reduction::Niggli)?;
    reduced.structure.p1()
}

fn standardized_primitive_cell(
    structure: &Structure,
    comparator: SpeciesComparator,
) -> Result<P1Cell, Error> {
    let conventional = standardize(structure, SymmetryOptions::default())?;
    let canonical = canonical_lattice(conventional.structure.lattice())?;
    let conventional = replace_lattice(&conventional.structure.p1()?, canonical)?;
    let primitive = primitive_cell(&conventional, SymmetryOptions::default())?;
    let reduced = reduce_structure(&primitive.structure, Reduction::Niggli)?;
    let canonical = canonical_lattice(reduced.structure.lattice())?;
    let cell = replace_lattice(&reduced.structure.p1()?, canonical)?;
    canonical_orientation(cell.p1()?, comparator)
}

/// Lattice with the same parameters in the standard orientation.
fn canonical_lattice(lattice: &Lattice) -> Result<Lattice, Error> {
    let parameters = lattice.parameters();
    Lattice::from_parameters(
        parameters[0],
        parameters[1],
        parameters[2],
        parameters[3],
        parameters[4],
        parameters[5],
    )
}

/// Rebuilds a P1 cell in another lattice, keeping fractional coordinates.
fn replace_lattice(cell: &P1Cell, lattice: Lattice) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new("standardized", lattice);
    for atom in &cell.atoms {
        builder = builder.site(
            Site::new(atom.label.clone(), &atom.element, atom.frac)?
                .with_occupancy(atom.occupancy)?,
        );
    }
    builder.build()
}

/// Picks the lexicographically first decoration among the lattice automorphisms.
fn canonical_orientation(cell: P1Cell, comparator: SpeciesComparator) -> Result<P1Cell, Error> {
    let automorphisms = lattice_automorphisms(&cell.lattice);
    let mut best: Option<CanonicalChoice> = None;
    for matrix in &automorphisms {
        let inverse = inverse_f64(*matrix)?;
        let mut decorated = Vec::with_capacity(cell.atoms.len());
        for atom in &cell.atoms {
            let frac = Lattice::wrap([
                atom.frac[0] * inverse[0][0]
                    + atom.frac[1] * inverse[1][0]
                    + atom.frac[2] * inverse[2][0],
                atom.frac[0] * inverse[0][1]
                    + atom.frac[1] * inverse[1][1]
                    + atom.frac[2] * inverse[2][1],
                atom.frac[0] * inverse[0][2]
                    + atom.frac[1] * inverse[1][2]
                    + atom.frac[2] * inverse[2][2],
            ]);
            decorated.push((atom, frac));
        }
        let mut signature: AtomSignature = decorated
            .iter()
            .map(|(atom, frac)| {
                (
                    signature_element(&atom.element, comparator),
                    quantize(frac[0]),
                    quantize(frac[1]),
                    quantize(frac[2]),
                )
            })
            .collect();
        signature.sort();
        if best
            .as_ref()
            .is_none_or(|(current, _)| signature < *current)
        {
            best = Some((signature, inverse));
        }
    }
    let Some((_, inverse)) = best else {
        return Ok(cell);
    };
    let mut builder = StructureBuilder::new("standardized", cell.lattice);
    for atom in &cell.atoms {
        let frac = Lattice::wrap([
            atom.frac[0] * inverse[0][0]
                + atom.frac[1] * inverse[1][0]
                + atom.frac[2] * inverse[2][0],
            atom.frac[0] * inverse[0][1]
                + atom.frac[1] * inverse[1][1]
                + atom.frac[2] * inverse[2][1],
            atom.frac[0] * inverse[0][2]
                + atom.frac[1] * inverse[1][2]
                + atom.frac[2] * inverse[2][2],
        ]);
        builder = builder.site(
            Site::new(atom.label.clone(), &atom.element, frac)?.with_occupancy(atom.occupancy)?,
        );
    }
    builder.build()?.p1()
}

fn quantize(value: f64) -> i64 {
    (value * 1e6).round() as i64
}

/// Species label used to choose between equivalent orientations.
fn signature_element(element: &str, comparator: SpeciesComparator) -> String {
    match comparator {
        SpeciesComparator::Exact => element.to_owned(),
        SpeciesComparator::Element => xtalforge_core::elements::element_by_symbol(element)
            .map(|info| info.symbol.clone())
            .unwrap_or_else(|_| element.to_owned()),
        SpeciesComparator::Ignore => String::new(),
    }
}

fn inverse_f64(matrix: [[i32; 3]; 3]) -> Result<[[f64; 3]; 3], Error> {
    let lattice = Lattice::new(matrix.map(|row| row.map(f64::from)))?;
    Ok(lattice.inverse())
}

/// Proper integer automorphisms of a lattice, identity first.
fn lattice_automorphisms(lattice: &Lattice) -> Vec<[[i32; 3]; 3]> {
    let base = lattice.matrix();
    let mut result = Vec::new();
    let mut entries = [0i32; 9];
    for code in 0..19683u32 {
        let mut value = code;
        for entry in &mut entries {
            *entry = (value % 3) as i32 - 1;
            value /= 3;
        }
        let matrix = [
            [entries[0], entries[1], entries[2]],
            [entries[3], entries[4], entries[5]],
            [entries[6], entries[7], entries[8]],
        ];
        if determinant(matrix) != 1 {
            continue;
        }
        if transformed_matches(matrix, base) {
            result.push(matrix);
        }
    }
    result
}

fn transformed_matches(matrix: [[i32; 3]; 3], base: [[f64; 3]; 3]) -> bool {
    let metric = gram(base);
    let scale = metric
        .iter()
        .flatten()
        .fold(1.0_f64, |acc, value| acc.max(value.abs()));
    for (row_index, row) in matrix.iter().enumerate() {
        for (column_index, column) in matrix.iter().enumerate() {
            let mut value = 0.0;
            for (first, first_entry) in row.iter().enumerate() {
                let Some(metric_row) = metric.get(first) else {
                    return false;
                };
                for (second, second_entry) in column.iter().enumerate() {
                    let base_value = metric_row.get(second).copied().unwrap_or(0.0);
                    value += f64::from(*first_entry) * base_value * f64::from(*second_entry);
                }
            }
            let target = metric
                .get(row_index)
                .and_then(|metric_row| metric_row.get(column_index))
                .copied()
                .unwrap_or(0.0);
            if (value - target).abs() / scale > AUTOMORPHISM_TOLERANCE {
                return false;
            }
        }
    }
    true
}

fn gram(matrix: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    std::array::from_fn(|row| {
        let first = matrix.get(row).copied().unwrap_or([0.0; 3]);
        std::array::from_fn(|column| {
            let second = matrix.get(column).copied().unwrap_or([0.0; 3]);
            first
                .iter()
                .zip(second)
                .map(|(first, second)| first * second)
                .sum()
        })
    })
}

fn determinant(matrix: [[i32; 3]; 3]) -> i32 {
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
}
