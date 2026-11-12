#![allow(clippy::indexing_slicing, clippy::needless_range_loop)] // fixed-size numeric kernels

//! Lattice matching and heterostructure stacking (Task 9.5.4).

use serde::{Deserialize, Serialize};

use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder};

/// One lattice match between two 2D cells.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct InterfaceMatch {
    /// Integer supercell of the film in its own 2D basis.
    pub film_matrix: [[i32; 2]; 2],
    /// Integer supercell of the substrate.
    pub substrate_matrix: [[i32; 2]; 2],
    /// Film in-plane vectors in Angstrom.
    pub film_vectors: [[f64; 2]; 2],
    /// Substrate in-plane vectors in Angstrom.
    pub substrate_vectors: [[f64; 2]; 2],
    /// Maximum relative length mismatch.
    pub strain: f64,
    /// Angle mismatch in degrees.
    pub angle_mismatch: f64,
    /// Area of the matched cell in square Angstrom.
    pub area: f64,
}

/// Options for [`lattice_matches`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InterfaceOptions {
    /// Maximum matrix entry.
    pub max_index: i32,
    /// Maximum area of the matched cell in square Angstrom.
    pub max_area: f64,
    /// Maximum relative length mismatch.
    pub length_tolerance: f64,
    /// Maximum angle mismatch in degrees.
    pub angle_tolerance: f64,
}

impl Default for InterfaceOptions {
    fn default() -> Self {
        Self {
            max_index: 3,
            max_area: 200.0,
            length_tolerance: 0.05,
            angle_tolerance: 5.0,
        }
    }
}

/// A supercell described by its integer matrix and two 2D vectors.
type Supercell = ([[i32; 2]; 2], [f64; 2], [f64; 2]);

/// Suggests commensurate 2D supercells for a film on a substrate.
///
/// The two in-plane lattice vectors are the first two rows of each lattice.
/// Matrices with entries up to `max_index` are enumerated, and matches are
/// ranked by strain and then area.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for invalid options or a degenerate
/// in-plane cell.
pub fn lattice_matches(
    film: &Structure,
    substrate: &Structure,
    options: &InterfaceOptions,
) -> Result<Vec<InterfaceMatch>, Error> {
    if options.max_index < 1 {
        return Err(Error::InvalidInput(
            "the interface matrix bound must be positive".to_owned(),
        ));
    }
    if !options.max_area.is_finite() || options.max_area <= 0.0 {
        return Err(Error::InvalidInput(
            "the interface area bound must be positive and finite".to_owned(),
        ));
    }
    let film_vectors = in_plane(film.lattice())?;
    let substrate_vectors = in_plane(substrate.lattice())?;
    let film_cells = supercells(film_vectors, options, options.max_area)?;
    let substrate_cells = supercells(substrate_vectors, options, options.max_area)?;
    let mut matches = Vec::new();
    for (film_matrix, film_u, film_v) in &film_cells {
        for (substrate_matrix, substrate_u, substrate_v) in &substrate_cells {
            let length_film = [norm2(*film_u), norm2(*film_v)];
            let length_substrate = [norm2(*substrate_u), norm2(*substrate_v)];
            if length_film.iter().any(|value| *value < 1e-6)
                || length_substrate.iter().any(|value| *value < 1e-6)
            {
                continue;
            }
            let strain = (0..2)
                .map(|index| {
                    ((length_film[index] - length_substrate[index]) / length_substrate[index]).abs()
                })
                .fold(0.0_f64, f64::max);
            if strain > options.length_tolerance {
                continue;
            }
            let angle_film = angle_between(*film_u, *film_v);
            let angle_substrate = angle_between(*substrate_u, *substrate_v);
            let angle_mismatch = (angle_film - angle_substrate).abs();
            if angle_mismatch > options.angle_tolerance {
                continue;
            }
            let area = length_film[0] * length_film[1] * angle_film.to_radians().sin().abs();
            matches.push(InterfaceMatch {
                film_matrix: *film_matrix,
                substrate_matrix: *substrate_matrix,
                film_vectors: [*film_u, *film_v],
                substrate_vectors: [*substrate_u, *substrate_v],
                strain,
                angle_mismatch,
                area,
            });
        }
    }
    matches.sort_by(|first, second| {
        first
            .strain
            .partial_cmp(&second.strain)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                first
                    .area
                    .partial_cmp(&second.area)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| {
                (first.film_matrix, first.substrate_matrix)
                    .cmp(&(second.film_matrix, second.substrate_matrix))
            })
    });
    if matches.len() > 64 {
        matches.truncate(64);
    }
    Ok(matches)
}

/// A stacked film/substrate heterostructure.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Heterostructure {
    /// Combined P1 structure.
    pub structure: Structure,
    /// In-plane strain applied to the film.
    pub strain: f64,
    /// Vacuum added above the film in Angstrom.
    pub vacuum: f64,
    /// Gap between the substrate and the film in Angstrom.
    pub gap: f64,
}

/// Stacks two slabs along the surface normal.
///
/// The substrate `a` and `b` define the combined cell; film atoms are mapped
/// onto that cell (which applies the in-plane strain) and placed `gap`
/// Angstrom above the highest substrate atom. `shift` is a fractional
/// translation of the film in the substrate plane.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty slab, a negative gap or
/// vacuum, or a builder error.
pub fn stack_slabs(
    film: &Structure,
    substrate: &Structure,
    gap: f64,
    shift: [f64; 3],
    vacuum: f64,
) -> Result<Heterostructure, Error> {
    if !gap.is_finite() || gap < 0.0 {
        return Err(Error::InvalidInput(
            "the interface gap must be non-negative and finite".to_owned(),
        ));
    }
    if !vacuum.is_finite() || vacuum < 0.0 {
        return Err(Error::InvalidInput(
            "the interface vacuum must be non-negative and finite".to_owned(),
        ));
    }
    if !shift.iter().all(|value| value.is_finite()) {
        return Err(Error::InvalidInput(
            "the interface shift must be finite".to_owned(),
        ));
    }
    let substrate_p1 = substrate.p1()?;
    let film_p1 = film.p1()?;
    if substrate_p1.atoms.is_empty() || film_p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "the interface needs two non-empty slabs".to_owned(),
        ));
    }
    let substrate_matrix = substrate_p1.lattice.matrix();
    let film_matrix = film_p1.lattice.matrix();
    let normal = substrate_p1.lattice.plane_normal([0.0, 0.0, 1.0])?;
    let substrate_cart: Vec<[f64; 3]> = substrate_p1
        .atoms
        .iter()
        .map(|atom| substrate_p1.lattice.to_cart(atom.frac))
        .collect();
    let film_cart: Vec<[f64; 3]> = film_p1
        .atoms
        .iter()
        .map(|atom| film_p1.lattice.to_cart(atom.frac))
        .collect();
    let substrate_height = substrate_cart
        .iter()
        .map(|point| dot(*point, normal))
        .fold(f64::NEG_INFINITY, f64::max)
        - substrate_cart
            .iter()
            .map(|point| dot(*point, normal))
            .fold(f64::INFINITY, f64::min);
    let film_min = film_cart
        .iter()
        .map(|point| dot(*point, normal))
        .fold(f64::INFINITY, f64::min);
    let film_max = film_cart
        .iter()
        .map(|point| dot(*point, normal))
        .fold(f64::NEG_INFINITY, f64::max);
    let film_height = film_max - film_min;
    let substrate_min = substrate_cart
        .iter()
        .map(|point| dot(*point, normal))
        .fold(f64::INFINITY, f64::min);
    let base = substrate_min;
    let total = substrate_height + gap + film_height + vacuum;
    // Build the combined cell with the substrate in-plane vectors and a
    // normal whose length is the total height.
    let in_plane_a = [
        substrate_matrix[0][0],
        substrate_matrix[0][1],
        substrate_matrix[0][2],
    ];
    let in_plane_b = [
        substrate_matrix[1][0],
        substrate_matrix[1][1],
        substrate_matrix[1][2],
    ];
    let c_vector = [normal[0] * total, normal[1] * total, normal[2] * total];
    let lattice = Lattice::new([in_plane_a, in_plane_b, c_vector])?;
    let inverse = lattice.inverse();
    let substrate_shift = [
        shift[0] * in_plane_a[0] + shift[1] * in_plane_b[0],
        shift[0] * in_plane_a[1] + shift[1] * in_plane_b[1],
        shift[0] * in_plane_a[2] + shift[1] * in_plane_b[2],
    ];
    let mut builder = StructureBuilder::new("interface", lattice);
    for atom in &substrate_p1.atoms {
        builder = builder.site(
            Site::new(atom.label.clone(), &atom.element, atom.frac)?
                .with_occupancy(atom.occupancy)?,
        );
    }
    let film_offset = substrate_height + gap;
    let mut strain = 0.0_f64;
    for (atom, cart) in film_p1.atoms.iter().zip(&film_cart) {
        let relative = [
            cart[0] + substrate_shift[0],
            cart[1] + substrate_shift[1],
            cart[2] - film_min + base + film_offset,
        ];
        let frac = [
            relative[0] * inverse[0][0] + relative[1] * inverse[0][1] + relative[2] * inverse[0][2],
            relative[0] * inverse[1][0] + relative[1] * inverse[1][1] + relative[2] * inverse[1][2],
            relative[0] * inverse[2][0] + relative[1] * inverse[2][1] + relative[2] * inverse[2][2],
        ];
        builder = builder.site(
            Site::new(atom.label.clone(), &atom.element, Lattice::wrap(frac))?
                .with_occupancy(atom.occupancy)?,
        );
    }
    let substrate_lengths = [norm3(in_plane_a), norm3(in_plane_b)];
    let film_lengths = [
        norm3([film_matrix[0][0], film_matrix[0][1], film_matrix[0][2]]),
        norm3([film_matrix[1][0], film_matrix[1][1], film_matrix[1][2]]),
    ];
    for (film_length, substrate_length) in film_lengths.iter().zip(substrate_lengths.iter()) {
        if *substrate_length > 1e-9 {
            strain = strain.max(((film_length - substrate_length) / substrate_length).abs());
        }
    }
    Ok(Heterostructure {
        structure: builder.build()?,
        strain,
        vacuum,
        gap,
    })
}

fn in_plane(lattice: &Lattice) -> Result<[[f64; 2]; 2], Error> {
    let matrix = lattice.matrix();
    let a = [matrix[0][0], matrix[0][1]];
    let b = [matrix[1][0], matrix[1][1]];
    if norm2(a) < 1e-9 || norm2(b) < 1e-9 {
        return Err(Error::InvalidInput(
            "the in-plane cell is degenerate".to_owned(),
        ));
    }
    Ok([a, b])
}

fn supercells(
    basis: [[f64; 2]; 2],
    options: &InterfaceOptions,
    max_area: f64,
) -> Result<Vec<Supercell>, Error> {
    let cell_area = (basis[0][0] * basis[1][1] - basis[0][1] * basis[1][0]).abs();
    if cell_area < 1e-9 {
        return Err(Error::InvalidInput(
            "the in-plane cell has zero area".to_owned(),
        ));
    }
    let limit = options.max_index;
    let mut result = Vec::new();
    for i in -limit..=limit {
        for j in -limit..=limit {
            let u = [
                f64::from(i) * basis[0][0] + f64::from(j) * basis[1][0],
                f64::from(i) * basis[0][1] + f64::from(j) * basis[1][1],
            ];
            if norm2(u) < 1e-6 {
                continue;
            }
            for k in -limit..=limit {
                for l in -limit..=limit {
                    let determinant = i * l - j * k;
                    if determinant <= 0 {
                        continue;
                    }
                    if f64::from(determinant) * cell_area > max_area {
                        continue;
                    }
                    let v = [
                        f64::from(k) * basis[0][0] + f64::from(l) * basis[1][0],
                        f64::from(k) * basis[0][1] + f64::from(l) * basis[1][1],
                    ];
                    if norm2(v) < 1e-6 {
                        continue;
                    }
                    result.push(([[i, j], [k, l]], u, v));
                }
            }
        }
    }
    result.sort_by(|first, second| {
        let area_first = (first.1[0] * first.2[1] - first.1[1] * first.2[0]).abs();
        let area_second = (second.1[0] * second.2[1] - second.1[1] * second.2[0]).abs();
        area_first
            .partial_cmp(&area_second)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| first.0.cmp(&second.0))
    });
    Ok(result)
}

fn angle_between(first: [f64; 2], second: [f64; 2]) -> f64 {
    let dot = first[0] * second[0] + first[1] * second[1];
    let cosine = dot / (norm2(first) * norm2(second));
    cosine.clamp(-1.0, 1.0).acos().to_degrees()
}

fn norm2(vector: [f64; 2]) -> f64 {
    (vector[0] * vector[0] + vector[1] * vector[1]).sqrt()
}

fn norm3(vector: [f64; 3]) -> f64 {
    (vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]).sqrt()
}

fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    first[0] * second[0] + first[1] * second[1] + first[2] * second[2]
}
