//! Deterministic 3D embedding of SMILES graphs (Tasks 6.1.2 and 6.1.3).
//!
//! This is a deliberately small distance-geometry embedder, not ETKDG. It
//! initializes the atoms on a seeded Fibonacci spiral, gives every bond an
//! element-pair ideal length, and relaxes the graph with a gradient descent on
//! a bond-length, bond-angle (1-3 distance) and non-bonded repulsion
//! objective. Ideal angles follow the hybridization implied by the bond
//! orders: 180, 120 or 109.47 degrees. Rings close because the
//! ring-closing bond is part of the objective; the result makes no
//! stereochemical claim and can be trapped in a local minimum for large or
//! highly fused graphs.
//!
//! Metal coordination is applied after the organic embedding: the metal's
//! neighbours are moved onto the vertices of an ideal polyhedron, and the rest
//! of the graph is left where the embedder put it. This is the same correction
//! Crystal-Builder applies to RDKit geometries (`build/coordination.py`), at a
//! reduced level of theory.

mod angles;
mod geometry;
mod hydrogens;
mod relax;

use xtalforge_core::elements::element_by_symbol;
use xtalforge_core::{Error, Structure};

use crate::smiles::SmilesMolecule;
use angles::angle_targets;
use geometry::{
    extra_direction, ideal_bond_length, ideal_bonds, is_metal, neighbours_of, radii,
    repulsion_pairs, shape_vertices, spiral,
};
use hydrogens::inferred_hydrogens;
use relax::relax;

/// Default value of [`EmbedOptions::bond_length`], Angstrom.
pub const DEFAULT_BOND_LENGTH: f64 = 1.5;

/// Padding around an embedded molecule, Angstrom.
const PADDING: f64 = 20.0;

/// Options for [`embed_smiles`] and [`apply_coordination`].
#[derive(Debug, Clone, PartialEq)]
pub struct EmbedOptions {
    /// Fallback bond length for element pairs without a table entry, Angstrom.
    pub bond_length: f64,
    /// Number of gradient steps.
    pub iterations: usize,
    /// Seed of the initial spiral and of its rotation.
    pub seed: u64,
}

impl Default for EmbedOptions {
    fn default() -> Self {
        Self {
            bond_length: DEFAULT_BOND_LENGTH,
            iterations: 1500,
            seed: 7,
        }
    }
}

/// Coordination polyhedron of a metal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CoordinationGeometry {
    /// Number of shape vertices.
    pub coordination: u8,
    /// Human-readable shape name.
    pub shape: &'static str,
}

/// Returns the usual coordination geometry of a metal, when tabulated.
///
/// The table is a small, documented set, not a replacement for a
/// crystallographic database: Fe, Mg, Ca, Al, Co, Mn and Cr are octahedral;
/// Zn is tetrahedral; Cu, Ni, Pt and Pd are square planar; Ag and Au are
/// linear.
#[must_use]
pub fn coordination_geometry(element: &str) -> Option<CoordinationGeometry> {
    let symbol = element_by_symbol(element).ok()?.symbol.as_str();
    let (coordination, shape) = match symbol {
        "Fe" | "Mg" | "Ca" | "Al" | "Co" | "Mn" | "Cr" => (6, "octahedral"),
        "Zn" => (4, "tetrahedral"),
        "Cu" | "Ni" | "Pt" | "Pd" => (4, "square planar"),
        "Ag" | "Au" => (2, "linear"),
        _ => return None,
    };
    Some(CoordinationGeometry {
        coordination,
        shape,
    })
}

/// Embeds a SMILES graph in three dimensions and returns a padded molecule.
///
/// The result is a P1 structure whose sites follow the SMILES atom order, so
/// bond `(a, b)` connects sites `a` and `b`. No explicit bonds are stored:
/// distance-based perception recovers them. The box is padded by 20 Angstrom on
/// every side.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for an empty graph, a non-positive or
/// non-finite bond length, zero iterations, an unknown element, a bond index
/// outside the graph, a non-positive bond order, or a valence that asks for
/// more hydrogens than the ideal geometry can place.
pub fn embed_smiles(graph: &SmilesMolecule, options: &EmbedOptions) -> Result<Structure, Error> {
    let points = embed_points(graph, options)?;
    molecule(graph, &points, false)
}

/// Embeds a SMILES graph and appends the hydrogens its bond orders imply.
///
/// Target valences are C 4, N 3, O 2, B 3, P 3 and S 2, shifted by the formal
/// charge. Aromatic bonds count as 1.5 and the used valence is rounded to
/// whole bonds, matching the hydrogen-filling rule of
/// [`crate::add_hydrogens`]. The heavy-atom embed is identical to
/// [`embed_smiles`]; the new H sites follow the graph atoms in atom order.
///
/// # Errors
///
/// Returns the errors of [`embed_smiles`] or a valence that asks for more
/// hydrogens than the ideal geometry can place.
pub fn embed_smiles_with_hydrogens(
    graph: &SmilesMolecule,
    options: &EmbedOptions,
) -> Result<Structure, Error> {
    let points = embed_points(graph, options)?;
    molecule(graph, &points, true)
}

/// Embeds a graph, then places the first metal atom's neighbours on `geometry`.
///
/// The metal itself keeps its embedded position. Each bonded neighbour is put
/// at one shape vertex, at the ideal metal-ligand distance, in SMILES bond
/// order. The rest of the graph is unchanged, so a ligand is not rigidly
/// moved with its donor atom.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for the errors of [`embed_smiles`], a graph
/// without a metal atom, a metal without neighbours, an unknown shape, a
/// coordination number that disagrees with the shape, or a metal geometry whose
/// coordination is zero.
///
/// # Warnings
///
/// A graph whose metal degree differs from the shape's coordination returns a
/// human-readable warning in the second tuple element. The first `degree`
/// neighbours still receive distinct vertices; extras use deterministic
/// fallback directions.
pub fn apply_coordination(
    graph: &SmilesMolecule,
    geometry: &CoordinationGeometry,
    options: &EmbedOptions,
) -> Result<(Structure, Vec<String>), Error> {
    validate_options(options)?;
    let metal = graph
        .atoms
        .iter()
        .position(|atom| is_metal(&atom.element))
        .ok_or_else(|| Error::InvalidInput("the SMILES graph has no metal atom".to_owned()))?;
    let vertices = shape_vertices(geometry)?;
    let neighbours = neighbours_of(graph, metal);
    if neighbours.is_empty() {
        return Err(Error::InvalidInput(
            "the metal atom has no bonded neighbours".to_owned(),
        ));
    }
    let mut warnings = Vec::new();
    if neighbours.len() != usize::from(geometry.coordination) {
        warnings.push(format!(
            "metal atom {metal} has {} neighbours but {} expects {}",
            neighbours.len(),
            geometry.shape,
            geometry.coordination
        ));
    }
    let mut points = embed_points(graph, options)?;
    let metal_point = points
        .get(metal)
        .copied()
        .ok_or_else(|| Error::InvalidInput("metal atom is outside the graph".to_owned()))?;
    let metal_element = graph
        .atoms
        .get(metal)
        .map(|atom| atom.element.as_str())
        .unwrap_or("C");
    for (index, (other, order)) in neighbours.iter().enumerate() {
        let direction = vertices
            .get(index)
            .copied()
            .unwrap_or_else(|| extra_direction(index));
        let Some(ligand) = graph.atoms.get(*other) else {
            continue;
        };
        let length = ideal_bond_length(metal_element, &ligand.element, *order, options);
        if let Some(slot) = points.get_mut(*other) {
            *slot = [
                metal_point[0] + direction[0] * length,
                metal_point[1] + direction[1] * length,
                metal_point[2] + direction[2] * length,
            ];
        }
    }
    Ok((molecule(graph, &points, false)?, warnings))
}

/// Returns the embedded Cartesian points, in graph atom order.
fn embed_points(graph: &SmilesMolecule, options: &EmbedOptions) -> Result<Vec<[f64; 3]>, Error> {
    validate_options(options)?;
    let count = graph.atoms.len();
    if count == 0 {
        return Err(Error::InvalidInput(
            "an embedding needs at least one atom".to_owned(),
        ));
    }
    let radii = radii(graph)?;
    let bonds = ideal_bonds(graph, options)?;
    let angles = angle_targets(graph, &bonds);
    let mut constrained = bonds.clone();
    constrained.extend(angles.iter().copied());
    let repulsion = repulsion_pairs(count, &constrained, &radii);
    let mut points = spiral(count, options);
    relax(&mut points, &bonds, &angles, &repulsion, options.iterations);
    Ok(points)
}

/// Wraps points as a `Structure::molecule` with the standard padding.
fn molecule(
    graph: &SmilesMolecule,
    points: &[[f64; 3]],
    infer_hydrogens: bool,
) -> Result<Structure, Error> {
    let mut named: Vec<(String, [f64; 3])> = graph
        .atoms
        .iter()
        .zip(points.iter())
        .map(|(atom, point)| (atom.element.clone(), *point))
        .collect();
    if infer_hydrogens {
        named.extend(inferred_hydrogens(graph, points)?);
    }
    Structure::molecule("molecule", &named, PADDING)
}

fn validate_options(options: &EmbedOptions) -> Result<(), Error> {
    if !options.bond_length.is_finite() || options.bond_length <= 0.0 {
        return Err(Error::InvalidInput(
            "embedding bond length must be positive and finite".to_owned(),
        ));
    }
    if options.iterations == 0 || options.iterations > 100_000 {
        return Err(Error::InvalidInput(format!(
            "embedding iterations must be in 1..=100000, got {}",
            options.iterations
        )));
    }
    Ok(())
}
