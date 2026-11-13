//! Structure-to-net contraction (Task 14.4.2).
//!
//! A crystal structure carries atoms; a net carries nodes and edges. This
//! module reduces one to the other in the three ways the rest of the library
//! needs:
//!
//! * [`NetNodeRule::Drawn`] reads the topology bonds stored on a P1 structure
//!   (the net the user drew over the chemistry).
//! * [`NetNodeRule::Atoms`] uses the perceived bond graph directly, optionally
//!   contracting every two-connected atom (the topological simplification).
//! * [`NetNodeRule::MetalClusters`] applies the "all-node" rule of the MOF
//!   deconstruction: a metal and its non-carbon coordination sphere, including
//!   hydrogens on those atoms, are one node; a carbon-containing unit bonded to
//!   exactly two nodes is one edge; a unit bonded to three or more nodes is
//!   itself a node. Terminal groups and free solvent are left out. The
//!   implementation lives in the private `contraction_metals` module.
//!
//! The output is a [`QuotientGraph`] with lattice shifts, so the existing
//! invariants, RCSR identification ([`crate::RcsrIndex`]) and vertex symbols
//! ([`mod@crate::vertex_symbols`]) apply unchanged.
//!
//! # Deliberate subset
//!
//! Node centroids are arithmetic means of unwrapped atom positions with the
//! first atom as the anchor. A cluster that wraps more than half a cell in two
//! directions at once is placed in the cell of that anchor, which is the same
//! convention [`xtalforge_io::cgd`] uses. Edge shifts are the nearest lattice
//! image of one centroid from the other, so the contraction is exact when two
//! bonded units are closer than half the cell along each axis.

#![allow(clippy::indexing_slicing)] // bounded local node and edge arrays

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{BondGraph, BondRules, Error, P1Cell, Site, Structure, StructureBuilder};

use crate::nets::{NetEdge, QuotientGraph};

use crate::contraction_simplify::strip_degree_two;

/// How a structure is reduced to net nodes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NetNodeRule {
    /// Use the topology bonds stored on a P1 structure.
    Drawn,
    /// Every bonded atom is a node.
    Atoms,
    /// Metals with their coordination spheres are nodes; linkers contract.
    #[default]
    MetalClusters,
}

/// Options for [`contract`].
#[derive(Debug, Clone, PartialEq)]
pub struct ContractionOptions {
    /// Node rule.
    pub rule: NetNodeRule,
    /// Contract two-connected nodes into one edge (topological simplification).
    /// Ignored by [`NetNodeRule::MetalClusters`], which already contracts
    /// two-connected linkers.
    pub contract_degree_two: bool,
    /// P1 expansion tolerance in Angstrom.
    pub tolerance: f64,
}

impl Default for ContractionOptions {
    fn default() -> Self {
        Self {
            rule: NetNodeRule::default(),
            contract_degree_two: false,
            tolerance: 0.05,
        }
    }
}

/// What a contracted node represents.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum NodeKind {
    /// One atom of an [`NetNodeRule::Atoms`] net.
    Atom,
    /// A metal with its coordination sphere.
    MetalCluster,
    /// A carbon-containing unit with three or more node connections.
    Linker,
    /// A node named by a drawn topology bond.
    Drawn,
}

/// One node of a contracted net.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ContractedNode {
    /// Atoms of the P1 structure collected into this node, sorted.
    pub atoms: Vec<usize>,
    /// Fractional centroid.
    pub centroid: [f64; 3],
    /// What the node represents.
    pub kind: NodeKind,
    /// Short label (`M`, `L`, an element symbol).
    pub label: String,
}

/// One edge of a contracted net.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct ContractedEdge {
    /// First node.
    pub a: usize,
    /// Second node.
    pub b: usize,
    /// Lattice shift applied to `b`.
    pub shift: [i32; 3],
    /// Linker atoms contracted into the edge, sorted; two atom indices for a
    /// direct bond and empty for a drawn topology edge.
    pub atoms: Vec<usize>,
}

/// A structure and the net it contracts to.
#[derive(Debug, Clone, PartialEq)]
pub struct ContractedNet {
    /// P1 structure the contraction was computed from.
    pub structure: Structure,
    /// Quotient graph; `nodes` and `edges` agree with the two fields below.
    pub net: QuotientGraph,
    /// Node descriptions in net order.
    pub nodes: Vec<ContractedNode>,
    /// Edge descriptions in net order.
    pub edges: Vec<ContractedEdge>,
    /// Node indices removed by two-connected contraction, in the numbering
    /// before the surviving nodes are renumbered.
    pub stripped: Vec<usize>,
}

impl ContractedNet {
    /// Coordination numbers and counts of the contracted net.
    #[must_use]
    pub fn invariants(&self) -> crate::NetInvariants {
        self.net.invariants()
    }

    /// Vertex symbols of the contracted net.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`crate::vertex_symbols::vertex_symbols`].
    pub fn vertex_symbols(
        &self,
        max_ring: usize,
    ) -> Result<Vec<crate::vertex_symbols::VertexSymbol>, Error> {
        crate::vertex_symbols::vertex_symbols(&self.net, max_ring)
    }

    /// Identifies the contracted net against the embedded RCSR index.
    ///
    /// # Errors
    ///
    /// Returns the errors of [`crate::RcsrIndex::identify`] or of loading the
    /// embedded index.
    pub fn identify(&self) -> Result<crate::rcsr::RcsrIdentification, Error> {
        crate::RcsrIndex::from_embedded()?.identify(&self.net)
    }
}

/// Contracts a structure to its underlying net.
///
/// # Errors
///
/// Returns [`Error::InvalidInput`] for a non-positive or non-finite tolerance,
/// a drawn topology on a non-P1 structure, a structure with no drawn topology
/// bonds, a structure with no metal for [`NetNodeRule::MetalClusters`], or a
/// contraction that leaves no node. Returns the errors of the P1 expansion,
/// bond perception or element lookup.
pub fn contract(
    structure: &Structure,
    options: &ContractionOptions,
) -> Result<ContractedNet, Error> {
    if !options.tolerance.is_finite() || options.tolerance <= 0.0 {
        return Err(Error::InvalidInput(
            "contraction tolerance must be positive and finite".to_owned(),
        ));
    }
    match options.rule {
        NetNodeRule::Drawn => {
            if !structure.is_p1() {
                return Err(Error::InvalidInput(
                    "a drawn topology needs a P1 structure; expand the structure first".to_owned(),
                ));
            }
            if structure.topology_bonds().is_empty() {
                return Err(Error::InvalidInput(
                    "the structure carries no topology bonds".to_owned(),
                ));
            }
            contract_drawn(structure, options)
        }
        NetNodeRule::Atoms => {
            let p1 = expansion(structure, options.tolerance)?;
            let bonds = xtalforge_core::perceive_bonds(&p1, BondRules::default())?;
            contract_atoms(&p1, &bonds, options)
        }
        NetNodeRule::MetalClusters => {
            let p1 = expansion(structure, options.tolerance)?;
            let bonds = xtalforge_core::perceive_bonds(&p1, BondRules::default())?;
            crate::contraction_metals::contract_metal_clusters(&p1, &bonds)
        }
    }
}

pub(crate) fn expansion(structure: &Structure, tolerance: f64) -> Result<Structure, Error> {
    let cell = P1Cell::expand(structure, tolerance)?;
    let mut builder = StructureBuilder::new(structure.name(), *structure.lattice()).hall("P 1");
    for atom in &cell.atoms {
        builder = builder.site(
            Site::new(atom.label.clone(), &atom.element, atom.frac)?
                .with_occupancy(atom.occupancy)?,
        );
    }
    builder.build()
}

fn contract_drawn(
    structure: &Structure,
    options: &ContractionOptions,
) -> Result<ContractedNet, Error> {
    let mut atoms = std::collections::BTreeSet::new();
    for bond in structure.topology_bonds() {
        atoms.insert(bond.i);
        atoms.insert(bond.j);
    }
    let picked: Vec<usize> = atoms.into_iter().collect();
    let position: BTreeMap<usize, usize> = picked
        .iter()
        .enumerate()
        .map(|(index, atom)| (*atom, index))
        .collect();
    let mut nodes = Vec::with_capacity(picked.len());
    let mut labels = Vec::with_capacity(picked.len());
    let mut descriptions = Vec::with_capacity(picked.len());
    for atom in &picked {
        let site = structure.sites().get(*atom).ok_or_else(|| {
            Error::InvalidInput("a topology bond references a missing atom".to_owned())
        })?;
        nodes.push(site.frac());
        labels.push(site.element().to_owned());
        descriptions.push(ContractedNode {
            atoms: vec![*atom],
            centroid: site.frac(),
            kind: NodeKind::Drawn,
            label: site.element().to_owned(),
        });
    }
    let mut edges = Vec::with_capacity(structure.topology_bonds().len());
    for bond in structure.topology_bonds() {
        let (Some(&a), Some(&b)) = (position.get(&bond.i), position.get(&bond.j)) else {
            return Err(Error::InvalidInput(
                "a topology bond references a missing atom".to_owned(),
            ));
        };
        edges.push(make_edge(a, b, bond.image, Vec::new()));
    }
    finish(
        structure.clone(),
        nodes,
        labels,
        edges,
        descriptions,
        options.contract_degree_two,
    )
}

fn contract_atoms(
    structure: &Structure,
    bonds: &BondGraph,
    options: &ContractionOptions,
) -> Result<ContractedNet, Error> {
    let n_atoms = structure.sites().len();
    let mut degree = vec![0usize; n_atoms];
    for bond in &bonds.bonds {
        if let Some(value) = degree.get_mut(bond.i) {
            *value += 1;
        }
        if bond.i != bond.j
            && let Some(value) = degree.get_mut(bond.j)
        {
            *value += 1;
        }
    }
    let active: Vec<usize> = (0..n_atoms).filter(|index| degree[*index] > 0).collect();
    if active.is_empty() {
        return Err(Error::InvalidInput(
            "the structure has no bonded atoms to contract".to_owned(),
        ));
    }
    let mut position = vec![usize::MAX; n_atoms];
    let mut nodes = Vec::with_capacity(active.len());
    let mut labels = Vec::with_capacity(active.len());
    let mut descriptions = Vec::with_capacity(active.len());
    for (vertex, atom) in active.iter().enumerate() {
        position[*atom] = vertex;
        let site = structure
            .sites()
            .get(*atom)
            .ok_or_else(|| Error::InvalidInput("a bond references a missing atom".to_owned()))?;
        nodes.push(site.frac());
        labels.push(site.element().to_owned());
        descriptions.push(ContractedNode {
            atoms: vec![*atom],
            centroid: site.frac(),
            kind: NodeKind::Atom,
            label: site.element().to_owned(),
        });
    }
    let mut edges = Vec::with_capacity(bonds.bonds.len());
    for bond in &bonds.bonds {
        let (Some(&first), Some(&second)) = (position.get(bond.i), position.get(bond.j)) else {
            continue;
        };
        if first == usize::MAX || second == usize::MAX {
            continue;
        }
        edges.push(make_edge(first, second, bond.image, vec![bond.i, bond.j]));
    }
    finish(
        structure.clone(),
        nodes,
        labels,
        edges,
        descriptions,
        options.contract_degree_two,
    )
}

pub(crate) fn finish(
    structure: Structure,
    mut nodes: Vec<[f64; 3]>,
    mut labels: Vec<String>,
    edges: Vec<ContractedEdge>,
    mut descriptions: Vec<ContractedNode>,
    strip: bool,
) -> Result<ContractedNet, Error> {
    let mut stripped = Vec::new();
    let mut edges = if strip {
        let (kept_nodes, kept_labels, kept_edges, removed) = strip_degree_two(nodes, labels, edges);
        nodes = kept_nodes;
        labels = kept_labels;
        stripped.extend(removed);
        kept_edges
    } else {
        edges
    };
    edges.sort();
    if nodes.is_empty() {
        return Err(Error::InvalidInput(
            "the structure contracts to no node".to_owned(),
        ));
    }
    if strip {
        descriptions = descriptions
            .into_iter()
            .enumerate()
            .filter(|(index, _)| !stripped.contains(index))
            .map(|(_, description)| description)
            .collect();
    }
    descriptions.truncate(nodes.len());
    for (index, description) in descriptions.iter_mut().enumerate() {
        if let Some(frac) = nodes.get(index) {
            description.centroid = *frac;
        }
    }
    let net_edges = edges
        .iter()
        .map(|edge| NetEdge {
            a: edge.a,
            b: edge.b,
            shift: edge.shift,
        })
        .collect();
    Ok(ContractedNet {
        structure,
        net: QuotientGraph {
            nodes,
            edges: net_edges,
            labels,
        },
        nodes: descriptions,
        edges,
        stripped,
    })
}

pub(crate) fn make_edge(a: usize, b: usize, shift: [i32; 3], atoms: Vec<usize>) -> ContractedEdge {
    let mut atoms = atoms;
    atoms.sort_unstable();
    atoms.dedup();
    let (a, b, shift) = if a <= b {
        (a, b, shift)
    } else {
        (b, a, shift.map(|value| -value))
    };
    let shift = if a == b {
        let reverse = shift.map(|value| -value);
        if reverse < shift { reverse } else { shift }
    } else {
        shift
    };
    ContractedEdge { a, b, shift, atoms }
}

pub(crate) fn lattice_shift(first: [f64; 3], second: [f64; 3]) -> [i32; 3] {
    [
        (first[0] - second[0]).round() as i32,
        (first[1] - second[1]).round() as i32,
        (first[2] - second[2]).round() as i32,
    ]
}

pub(crate) fn centroid(structure: &Structure, atoms: &[usize]) -> Result<[f64; 3], Error> {
    let anchor = atoms
        .first()
        .and_then(|atom| structure.sites().get(*atom))
        .map(Site::frac)
        .ok_or_else(|| Error::InvalidInput("a node has no atom".to_owned()))?;
    let mut sum = [0.0; 3];
    for &atom in atoms {
        let site = structure
            .sites()
            .get(atom)
            .ok_or_else(|| Error::InvalidInput("a node references a missing atom".to_owned()))?;
        for axis in 0..3 {
            let delta = site.frac()[axis] - anchor[axis];
            sum[axis] += anchor[axis] + delta - delta.round();
        }
    }
    let count = atoms.len() as f64;
    let mut out = [
        (sum[0] / count).rem_euclid(1.0),
        (sum[1] / count).rem_euclid(1.0),
        (sum[2] / count).rem_euclid(1.0),
    ];
    for value in &mut out {
        if (*value - 1.0).abs() < 1e-12 {
            *value = 0.0;
        }
    }
    Ok(out)
}

pub(crate) fn adjacency(bonds: &BondGraph, n_atoms: usize) -> Vec<Vec<usize>> {
    let mut out = vec![Vec::new(); n_atoms];
    for bond in &bonds.bonds {
        if let Some(list) = out.get_mut(bond.i) {
            list.push(bond.j);
        }
        if let Some(list) = out.get_mut(bond.j)
            && bond.i != bond.j
        {
            list.push(bond.i);
        }
    }
    out
}

pub(crate) fn union(parent: &mut [usize], first: usize, second: usize) {
    let first_root = find_root(parent, first);
    let second_root = find_root(parent, second);
    if first_root != second_root
        && let Some(slot) = parent.get_mut(second_root)
    {
        *slot = first_root;
    }
}

pub(crate) fn find_root(parent: &mut [usize], mut node: usize) -> usize {
    while parent.get(node).copied().unwrap_or(node) != node {
        let next = parent.get(node).copied().unwrap_or(node);
        let grandparent = parent.get(next).copied().unwrap_or(next);
        if let Some(slot) = parent.get_mut(node) {
            *slot = grandparent;
        }
        node = next;
    }
    node
}
