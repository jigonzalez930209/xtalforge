#![allow(clippy::indexing_slicing)] // bounded indices from validated loops and fixed-size vectors

//! Periodic crystal graphs (Task 21.1.1).
//!
//! The graph uses the P1 cell in a canonical order (atomic number, then
//! fractional position, then label). Every edge stores the lattice image of
//! the target and the Cartesian vector from source to target, so a GNN can
//! reconstruct the geometry without the cell.
//!
//! The neighbour list comes from the exact linked-cell search of Task 1.10.
//! Both directions of every pair are included; `max_neighbors` keeps the
//! shortest outgoing edges per node and marks the graph as truncated.
//!
//! Arrays are plain `Vec`s in Rust and `Vec<Vec<f64>>`/`Vec<Vec<usize>>` in
//! Python and WASM, which the adapters convert to NumPy or typed arrays.

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Lattice, Site, Structure, StructureBuilder, neighbors};

/// Options for [`crystal_graph`].
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GraphOptions {
    /// Neighbour cutoff in Angstrom.
    pub cutoff: f64,
    /// Keep at most this many outgoing edges per node; `None` keeps all.
    pub max_neighbors: Option<usize>,
}

impl Default for GraphOptions {
    fn default() -> Self {
        Self {
            cutoff: 5.0,
            max_neighbors: None,
        }
    }
}

type Edge = ((usize, usize, [i32; 3]), [f64; 3], f64);

/// A directed periodic graph.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CrystalGraph {
    /// Element symbol for each node.
    pub elements: Vec<String>,
    /// Atomic number for each node.
    pub atomic_numbers: Vec<u8>,
    /// Fractional position of each node in the P1 cell.
    pub positions: Vec<[f64; 3]>,
    /// Lattice matrix, rows as vectors.
    pub lattice: [[f64; 3]; 3],
    /// Source node of each edge.
    pub edge_from: Vec<usize>,
    /// Target node of each edge.
    pub edge_to: Vec<usize>,
    /// Lattice image of the target.
    pub edge_image: Vec<[i32; 3]>,
    /// Cartesian vector from source to target in Angstrom.
    pub edge_vectors: Vec<[f64; 3]>,
    /// Distance in Angstrom for each edge.
    pub distances: Vec<f64>,
    /// Whether `max_neighbors` removed any edge.
    pub truncated: bool,
}

impl CrystalGraph {
    /// Number of nodes.
    #[must_use]
    pub fn num_nodes(&self) -> usize {
        self.elements.len()
    }

    /// Number of edges.
    #[must_use]
    pub fn num_edges(&self) -> usize {
        self.edge_from.len()
    }

    /// Edge pairs as `[source, target]`.
    #[must_use]
    pub fn edge_index(&self) -> Vec<[usize; 2]> {
        self.edge_from
            .iter()
            .zip(&self.edge_to)
            .map(|(from, to)| [*from, *to])
            .collect()
    }

    /// Serializes the graph to JSON.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when serialization fails.
    pub fn to_json(&self) -> Result<String, Error> {
        serde_json::to_string(self)
            .map_err(|error| Error::InvalidInput(format!("could not serialize graph: {error}")))
    }
}

/// Builds the periodic graph of a structure.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-positive cutoff or an empty
/// structure, and the neighbour or builder errors.
pub fn crystal_graph(structure: &Structure, options: &GraphOptions) -> Result<CrystalGraph, Error> {
    if !options.cutoff.is_finite() || options.cutoff <= 0.0 {
        return Err(Error::InvalidInput(
            "the graph cutoff must be positive and finite".to_owned(),
        ));
    }
    if options.max_neighbors == Some(0) {
        return Err(Error::InvalidInput(
            "max_neighbors must be at least one".to_owned(),
        ));
    }
    let canonical = canonical_p1(structure)?;
    let p1 = canonical.p1()?;
    if p1.atoms.is_empty() {
        return Err(Error::InvalidInput(
            "a graph needs at least one atom".to_owned(),
        ));
    }
    let raw = neighbors(&canonical, options.cutoff)?;
    let mut edges: Vec<Edge> = Vec::new();
    for neighbor in raw {
        let from = p1
            .atoms
            .get(neighbor.i)
            .map(|atom| atom.frac)
            .unwrap_or([0.0; 3]);
        let to = p1
            .atoms
            .get(neighbor.j)
            .map(|atom| atom.frac)
            .unwrap_or([0.0; 3]);
        let image = neighbor.image;
        let target = [
            to[0] + image[0] as f64,
            to[1] + image[1] as f64,
            to[2] + image[2] as f64,
        ];
        let vector = p1.lattice.to_cart([
            target[0] - from[0],
            target[1] - from[1],
            target[2] - from[2],
        ]);
        edges.push(((neighbor.i, neighbor.j, image), vector, neighbor.distance));
        let reverse_image = [-image[0], -image[1], -image[2]];
        edges.push((
            (neighbor.j, neighbor.i, reverse_image),
            [-vector[0], -vector[1], -vector[2]],
            neighbor.distance,
        ));
    }
    edges.sort_by(|left, right| {
        left.0
            .0
            .cmp(&right.0.0)
            .then_with(|| left.0.1.cmp(&right.0.1))
            .then_with(|| left.0.2.cmp(&right.0.2))
    });
    let mut truncated = false;
    if let Some(limit) = options.max_neighbors {
        let mut counts = vec![0usize; p1.atoms.len()];
        let mut kept = Vec::with_capacity(edges.len());
        for edge in edges {
            let source = edge.0.0;
            if counts.get(source).copied().unwrap_or(0) < limit {
                counts[source] += 1;
                kept.push(edge);
            } else {
                truncated = true;
            }
        }
        edges = kept;
    }
    let elements: Vec<String> = p1.atoms.iter().map(|atom| atom.element.clone()).collect();
    let mut atomic_numbers = Vec::with_capacity(elements.len());
    for element in &elements {
        atomic_numbers.push(xtalforge_core::elements::element_by_symbol(element)?.z);
    }
    let positions: Vec<[f64; 3]> = p1.atoms.iter().map(|atom| atom.frac).collect();
    let mut edge_from = Vec::with_capacity(edges.len());
    let mut edge_to = Vec::with_capacity(edges.len());
    let mut edge_image = Vec::with_capacity(edges.len());
    let mut edge_vectors = Vec::with_capacity(edges.len());
    let mut distances = Vec::with_capacity(edges.len());
    for ((from, to, image), vector, distance) in edges {
        edge_from.push(from);
        edge_to.push(to);
        edge_image.push(image);
        edge_vectors.push(vector);
        distances.push(distance);
    }
    Ok(CrystalGraph {
        elements,
        atomic_numbers,
        positions,
        lattice: p1.lattice.matrix(),
        edge_from,
        edge_to,
        edge_image,
        edge_vectors,
        distances,
        truncated,
    })
}

/// Builds a canonical-order P1 structure.
///
/// # Errors
///
/// Returns the P1 or builder errors.
pub fn canonical_p1(structure: &Structure) -> Result<Structure, Error> {
    let p1 = structure.p1()?;
    let mut atoms: Vec<&xtalforge_core::P1Atom> = p1.atoms.iter().collect();
    atoms.sort_by(|left, right| {
        let left_z = xtalforge_core::elements::element_by_symbol(&left.element)
            .map(|element| element.z)
            .unwrap_or(0);
        let right_z = xtalforge_core::elements::element_by_symbol(&right.element)
            .map(|element| element.z)
            .unwrap_or(0);
        left_z
            .cmp(&right_z)
            .then_with(|| left.frac[0].total_cmp(&right.frac[0]))
            .then_with(|| left.frac[1].total_cmp(&right.frac[1]))
            .then_with(|| left.frac[2].total_cmp(&right.frac[2]))
            .then_with(|| left.label.cmp(&right.label))
    });
    let mut builder = StructureBuilder::new(structure.name(), p1.lattice);
    for atom in atoms {
        builder = builder.site(Site::new(&atom.label, &atom.element, atom.frac)?);
    }
    builder.build()
}

/// The graph as Python-style arrays.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GraphArrays {
    /// Node atomic numbers.
    pub atomic_numbers: Vec<u8>,
    /// Node fractional positions.
    pub positions: Vec<[f64; 3]>,
    /// `[source, target]` pairs.
    pub edge_index: Vec<[usize; 2]>,
    /// Edge vectors.
    pub edge_vectors: Vec<[f64; 3]>,
    /// Edge distances.
    pub distances: Vec<f64>,
    /// Edge lattice images.
    pub edge_images: Vec<[i32; 3]>,
    /// Lattice matrix.
    pub lattice: [[f64; 3]; 3],
}

impl CrystalGraph {
    /// Arrays ready for NumPy or typed-array conversion.
    #[must_use]
    pub fn arrays(&self) -> GraphArrays {
        GraphArrays {
            atomic_numbers: self.atomic_numbers.clone(),
            positions: self.positions.clone(),
            edge_index: self.edge_index(),
            edge_vectors: self.edge_vectors.clone(),
            distances: self.distances.clone(),
            edge_images: self.edge_image.clone(),
            lattice: self.lattice,
        }
    }
}

/// Builds a graph and returns its arrays.
///
/// # Errors
///
/// Returns the errors of [`crystal_graph`].
pub fn crystal_graph_arrays(
    structure: &Structure,
    options: &GraphOptions,
) -> Result<GraphArrays, Error> {
    Ok(crystal_graph(structure, options)?.arrays())
}

/// Builds a graph from an explicit lattice and fractional positions.
///
/// # Errors
///
/// Returns the builder errors.
pub fn graph_from_positions(
    lattice: Lattice,
    elements: &[String],
    positions: &[[f64; 3]],
    options: &GraphOptions,
) -> Result<CrystalGraph, Error> {
    if elements.len() != positions.len() {
        return Err(Error::InvalidInput(
            "elements and positions must have the same length".to_owned(),
        ));
    }
    let mut builder = StructureBuilder::new("graph", lattice);
    for (index, (element, position)) in elements.iter().zip(positions).enumerate() {
        builder = builder.site(Site::new(
            format!("{element}{}", index + 1),
            element,
            *position,
        )?);
    }
    crystal_graph(&builder.build()?, options)
}
