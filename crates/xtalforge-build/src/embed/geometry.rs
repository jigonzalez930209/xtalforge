//! Geometry helpers for the reduced distance-geometry embedder.

use std::collections::BTreeSet;

use xtalforge_core::Error;
use xtalforge_core::elements::element_by_symbol;

use super::EmbedOptions;
use crate::smiles::SmilesMolecule;

/// Golden angle, radians.
pub(super) const GOLDEN_ANGLE: f64 = 2.399_963_229_728_653;
/// Below this separation two atoms have no usable direction.
pub(super) const DEGENERATE: f64 = 1e-9;
/// Fraction of the covalent-radius sum used as the repulsion distance.
const REPULSION_SCALE: f64 = 0.80;

/// Covalent radii for every graph atom, in atom order.
pub(super) fn radii(graph: &SmilesMolecule) -> Result<Vec<f64>, Error> {
    graph
        .atoms
        .iter()
        .map(|atom| Ok(element_by_symbol(&atom.element)?.covalent_radius))
        .collect()
}

/// Ideal bond lengths for every graph bond.
pub(super) fn ideal_bonds(
    graph: &SmilesMolecule,
    options: &EmbedOptions,
) -> Result<Vec<(usize, usize, f64)>, Error> {
    let mut bonds = Vec::with_capacity(graph.bonds.len());
    for bond in &graph.bonds {
        let first = graph
            .atoms
            .get(bond.a)
            .ok_or_else(|| Error::InvalidInput(format!("bond atom {} is missing", bond.a)))?;
        let second = graph
            .atoms
            .get(bond.b)
            .ok_or_else(|| Error::InvalidInput(format!("bond atom {} is missing", bond.b)))?;
        if !bond.order.is_finite() || bond.order <= 0.0 {
            return Err(Error::InvalidInput(
                "bond order must be positive and finite".to_owned(),
            ));
        }
        bonds.push((
            bond.a,
            bond.b,
            ideal_bond_length(&first.element, &second.element, bond.order, options),
        ));
    }
    Ok(bonds)
}

/// Ideal length of one bond, Angstrom.
pub(super) fn ideal_bond_length(
    first: &str,
    second: &str,
    order: f64,
    options: &EmbedOptions,
) -> f64 {
    let base = pair_length(first, second).unwrap_or(options.bond_length);
    let factor = if (order - 1.5).abs() < 1e-6 {
        0.93
    } else if (order - 2.0).abs() < 1e-6 {
        0.87
    } else if (order - 3.0).abs() < 1e-6 {
        0.78
    } else {
        1.0
    };
    base * factor
}

/// Single-bond reference length for an element pair, Angstrom.
fn pair_length(first: &str, second: &str) -> Option<f64> {
    let mut pair = [first, second];
    pair.sort_unstable();
    match pair {
        ["C", "C"] => Some(1.54),
        ["C", "N"] => Some(1.47),
        ["C", "O"] => Some(1.43),
        ["C", "H"] => Some(1.09),
        ["C", "F"] => Some(1.35),
        ["C", "Cl"] => Some(1.77),
        ["C", "Br"] => Some(1.94),
        ["C", "I"] => Some(2.14),
        ["C", "S"] => Some(1.82),
        ["N", "N"] => Some(1.45),
        ["N", "O"] => Some(1.40),
        ["O", "O"] => Some(1.48),
        ["H", "H"] => Some(0.74),
        ["H", "N"] => Some(1.01),
        ["O", "H"] => Some(0.96),
        ["N", "S"] => Some(1.68),
        ["O", "S"] => Some(1.58),
        _ => None,
    }
}

/// Initial deterministic spiral, one point per atom.
pub(super) fn spiral(count: usize, options: &EmbedOptions) -> Vec<[f64; 3]> {
    let phase = (options.seed % 360) as f64 * std::f64::consts::PI / 180.0;
    let scale = options.bond_length * (count as f64).cbrt().max(0.5);
    (0..count)
        .map(|index| {
            let fraction = (index as f64 + 0.5) / count as f64;
            let z = 1.0 - 2.0 * fraction;
            let ring = (1.0 - z * z).max(0.0).sqrt();
            let angle = GOLDEN_ANGLE * index as f64 + phase;
            [
                scale * ring * angle.cos(),
                scale * ring * angle.sin(),
                scale * z,
            ]
        })
        .collect()
}

/// Pairs that are not bonded and the distance below which they repel.
pub(super) fn repulsion_pairs(
    count: usize,
    bonds: &[(usize, usize, f64)],
    radii: &[f64],
) -> Vec<(usize, usize, f64)> {
    let bonded: BTreeSet<(usize, usize)> = bonds
        .iter()
        .map(|&(first, second, _)| {
            if first <= second {
                (first, second)
            } else {
                (second, first)
            }
        })
        .collect();
    let mut pairs = Vec::new();
    for first in 0..count {
        for second in (first + 1)..count {
            if bonded.contains(&(first, second)) {
                continue;
            }
            let radius = radii.get(first).copied().unwrap_or(0.76)
                + radii.get(second).copied().unwrap_or(0.76);
            pairs.push((first, second, REPULSION_SCALE * radius));
        }
    }
    pairs
}

/// Unit direction of a pair, with a deterministic fallback for coincident atoms.
pub(super) fn unit(delta: [f64; 3], first: usize, second: usize) -> [f64; 3] {
    let [dx, dy, dz] = delta;
    let distance = (dx * dx + dy * dy + dz * dz).sqrt();
    if distance > DEGENERATE {
        [dx / distance, dy / distance, dz / distance]
    } else {
        let angle = (first as f64 + 1.0) * GOLDEN_ANGLE + second as f64;
        [angle.cos(), angle.sin(), 0.31]
    }
}

/// Whether an element is a metal in the reference table.
pub(super) fn is_metal(element: &str) -> bool {
    element_by_symbol(element).is_ok_and(|entry| entry.is_metal)
}

/// Bonded neighbours of one atom, in graph bond order.
pub(super) fn neighbours_of(graph: &SmilesMolecule, index: usize) -> Vec<(usize, f64)> {
    let mut neighbours = Vec::new();
    for bond in &graph.bonds {
        if bond.a == index {
            neighbours.push((bond.b, bond.order));
        } else if bond.b == index {
            neighbours.push((bond.a, bond.order));
        }
    }
    neighbours
}

/// Vertices of an ideal coordination polyhedron as unit vectors.
pub(super) fn shape_vertices(
    geometry: &super::CoordinationGeometry,
) -> Result<Vec<[f64; 3]>, Error> {
    let vertices: Vec<[f64; 3]> = match geometry.shape {
        "linear" => vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
        "trigonal planar" => vec![
            [1.0, 0.0, 0.0],
            [-0.5, 0.866_025_403_784_438_6, 0.0],
            [-0.5, -0.866_025_403_784_438_6, 0.0],
        ],
        "tetrahedral" => {
            let value = 1.0 / 3.0_f64.sqrt();
            vec![
                [value, value, value],
                [value, -value, -value],
                [-value, value, -value],
                [-value, -value, value],
            ]
        }
        "square planar" => vec![
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
        ],
        "trigonal bipyramidal" => vec![
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
            [1.0, 0.0, 0.0],
            [-0.5, 0.866_025_403_784_438_6, 0.0],
            [-0.5, -0.866_025_403_784_438_6, 0.0],
        ],
        "square pyramidal" => vec![
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
        ],
        "octahedral" => vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ],
        other => {
            return Err(Error::InvalidInput(format!(
                "unknown coordination shape {other:?}"
            )));
        }
    };
    if vertices.len() != usize::from(geometry.coordination) {
        return Err(Error::InvalidInput(format!(
            "coordination {} does not match shape {:?} with {} vertices",
            geometry.coordination,
            geometry.shape,
            vertices.len()
        )));
    }
    Ok(vertices)
}

/// Deterministic direction for a neighbour beyond the polyhedron.
pub(super) fn extra_direction(index: usize) -> [f64; 3] {
    let angle = GOLDEN_ANGLE * (index as f64 + 1.0);
    [angle.cos(), angle.sin(), 0.5]
}
