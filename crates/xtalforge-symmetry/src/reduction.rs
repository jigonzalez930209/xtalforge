//! Lattice reduction and structure-level reduction (Task 2.6).
//!
//! moyo 0.22 exposes `Lattice::niggli_reduce`, `delaunay_reduce` and
//! `minkowski_reduce` publicly, so no vendored port is needed.

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

/// Reduction algorithm.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Reduction {
    /// Niggli reduction.
    Niggli,
    /// Delaunay reduction.
    Delaunay,
    /// Minkowski reduction.
    Minkowski,
}

/// A reduced structure and the row transformation that produced it.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReducedStructure {
    /// P1 structure in the reduced cell.
    pub structure: Structure,
    /// Row transformation: `new_lattice = P * old_lattice`.
    pub transformation: [[f64; 3]; 3],
    /// Determinant of the transformation (`+1` or `-1`).
    pub determinant: f64,
}

/// Reduces a lattice with the chosen algorithm.
///
/// # Errors
///
/// Returns [`Error::Symmetry`] when moyo refuses the reduction.
pub fn reduce_lattice(
    lattice: &Lattice,
    reduction: Reduction,
) -> Result<(Lattice, [[f64; 3]; 3]), Error> {
    let moyo_lattice = moyo::base::Lattice::from_basis(lattice.matrix());
    let (reduced, transform) = match reduction {
        Reduction::Niggli => moyo_lattice.niggli_reduce(),
        Reduction::Delaunay => moyo_lattice.delaunay_reduce(),
        Reduction::Minkowski => moyo_lattice.minkowski_reduce(),
    }
    .map_err(|error| Error::Symmetry(format!("lattice reduction failed: {error}")))?;
    let rows: [[f64; 3]; 3] = reduced.basis.transpose().into();
    let transform: [[f64; 3]; 3] = transform.map(f64::from).into();
    // basis' = basis * T, so rows' = T^T * rows.
    let [[a, b, c], [d, e, f], [g, h, i]] = transform;
    let row_transform = [[a, d, g], [b, e, h], [c, f, i]];
    let reduced = Lattice::new(rows)?;
    Ok((reduced, row_transform))
}

/// Reduces a structure to a P1 cell in the reduced lattice.
///
/// # Errors
///
/// Returns the errors of [`reduce_lattice`] or a builder error.
pub fn reduce_structure(
    structure: &Structure,
    reduction: Reduction,
) -> Result<ReducedStructure, Error> {
    let (reduced_lattice, transformation) = reduce_lattice(structure.lattice(), reduction)?;
    let determinant = determinant(transformation);
    let p1 = structure.p1()?;
    let mut sites: Vec<Site> = Vec::new();
    for atom in &p1.atoms {
        let frac = Lattice::wrap(reduced_lattice.to_frac(p1.lattice.to_cart(atom.frac)));
        let duplicate = sites.iter().any(|existing| {
            existing.element() == atom.element
                && reduced_lattice.distance(existing.frac(), frac) <= 1e-6
        });
        if !duplicate {
            sites.push(
                Site::new(atom.label.clone(), &atom.element, frac)?
                    .with_occupancy(atom.occupancy)?,
            );
        }
    }
    let mut builder = StructureBuilder::new(structure.name(), reduced_lattice);
    for site in sites {
        builder = builder.site(site);
    }
    for (key, value) in structure.metadata() {
        builder = builder.metadata(key.clone(), value.clone());
    }
    Ok(ReducedStructure {
        structure: builder.build()?,
        transformation,
        determinant,
    })
}

/// Whether two lattices are equivalent after Niggli reduction.
///
/// # Errors
///
/// Returns the errors of [`reduce_lattice`].
pub fn lattices_equivalent(
    first: &Lattice,
    second: &Lattice,
    tolerance: f64,
) -> Result<bool, Error> {
    let (first, _) = reduce_lattice(first, Reduction::Niggli)?;
    let (second, _) = reduce_lattice(second, Reduction::Niggli)?;
    let mut first_parameters = first.parameters();
    let mut second_parameters = second.parameters();
    // Compare lengths and angles as multisets: Niggli cells may be reported
    // with a different row order.
    first_parameters[..3].sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    second_parameters[..3].sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    first_parameters[3..].sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    second_parameters[3..].sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    Ok(first_parameters
        .iter()
        .zip(second_parameters.iter())
        .all(|(first, second)| (first - second).abs() <= tolerance))
}

fn determinant(matrix: [[f64; 3]; 3]) -> f64 {
    let [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
}
