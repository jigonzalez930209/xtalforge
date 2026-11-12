#![allow(clippy::indexing_slicing)] // fixed 3x3 and bounded graph arrays

//! Deconstruction of a framework into nodes and linkers (Task 10.2.1).

use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};
use xtalforge_analysis::nets::{NetEdge, QuotientGraph};
use xtalforge_core::{Error, Structure, element_info};

use crate::platform::util::{bond_adjacency, components, element_counts, formula, p1_and_bonds};

/// Kind of a building unit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum UnitKind {
    /// Metal cluster with its coordinating atoms.
    Node,
    /// Organic link between nodes.
    Linker,
    /// Free solvent or counter-ion fragment.
    Solvent,
}

/// One connection from a unit to another unit.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct UnitConnection {
    /// Parent P1 atom index on this unit.
    pub atom: usize,
    /// Parent P1 atom index on the other unit.
    pub other: usize,
    /// Index of the other unit in [`Deconstruction::units`].
    pub other_unit: usize,
    /// Lattice image of the other atom.
    pub image: [i32; 3],
}

/// One deconstructed building unit.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BuildingUnit {
    /// Node, linker or solvent.
    pub kind: UnitKind,
    /// Parent P1 atom indices, sorted.
    pub atoms: Vec<usize>,
    /// Chemical formula.
    pub formula: String,
    /// Element counts.
    pub elements: BTreeMap<String, usize>,
    /// Connections to other units.
    pub connections: Vec<UnitConnection>,
    /// Centroid in fractional coordinates.
    pub centroid: [f64; 3],
    /// Longest connection-to-connection distance in Angstrom.
    pub span: f64,
    /// Metal elements of a node.
    pub metals: Vec<String>,
}

/// The deconstructed framework and its underlying net.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Deconstruction {
    /// P1 structure used for the analysis.
    pub structure: Structure,
    /// Nodes, linkers and solvent fragments.
    pub units: Vec<BuildingUnit>,
    /// Contracted periodic net.
    pub net: QuotientGraph,
    /// RCSR names with the same invariants.
    pub topology: Vec<String>,
    /// Nearest RCSR names when no exact match exists.
    pub nearest_topologies: Vec<String>,
    /// Number of interpenetrating framework components.
    pub interpenetration: usize,
}

/// Deconstructs a framework into metal nodes and organic linkers.
///
/// The "all-node" rule assigns every metal atom and every non-carbon atom
/// directly bonded to a metal (plus its hydrogens) to a node. Metals are
/// joined into one node when they share a non-carbon bridge. The remaining
/// main-framework atoms form linkers. A two-connected linker becomes an edge
/// of the contracted net; higher-connected units become vertices.
///
/// # Errors
///
/// Returns the errors of the P1 expansion, bond perception, element lookup or
/// the RCSR index.
pub fn deconstruct(structure: &Structure) -> Result<Deconstruction, Error> {
    let (p1, graph) = p1_and_bonds(structure)?;
    let adjacency = bond_adjacency(&graph, p1.sites().len());
    let groups = components(&adjacency);
    let main = groups
        .iter()
        .enumerate()
        .max_by_key(|(_, group)| group.len())
        .map(|(index, _)| index)
        .unwrap_or(0);
    let mut node_atoms: BTreeSet<usize> = BTreeSet::new();
    let mut metals: Vec<usize> = Vec::new();
    for (index, site) in p1.sites().iter().enumerate() {
        let info = element_info(site.element())?;
        if info.is_metal {
            metals.push(index);
            node_atoms.insert(index);
        }
    }
    let metal_set: BTreeSet<usize> = metals.iter().copied().collect();
    for &metal in &metals {
        for &neighbour in adjacency.get(metal).map(Vec::as_slice).unwrap_or(&[]) {
            let Some(site) = p1.sites().get(neighbour) else {
                continue;
            };
            if site.element() != "C" {
                node_atoms.insert(neighbour);
                for &second in adjacency.get(neighbour).map(Vec::as_slice).unwrap_or(&[]) {
                    if p1
                        .sites()
                        .get(second)
                        .is_some_and(|site| site.element() == "H")
                    {
                        node_atoms.insert(second);
                    }
                }
            }
        }
    }
    // Join metals that share a non-carbon bridge.
    let mut metal_parent: Vec<usize> = (0..p1.sites().len()).collect();
    for &bridge in &node_atoms {
        if metal_set.contains(&bridge) {
            continue;
        }
        let bonded: Vec<usize> = adjacency
            .get(bridge)
            .map(Vec::as_slice)
            .unwrap_or(&[])
            .iter()
            .copied()
            .filter(|index| metal_set.contains(index))
            .collect();
        if let Some(&first) = bonded.first() {
            for &other in bonded.iter().skip(1) {
                let first_root = find_root(&mut metal_parent, first);
                let other_root = find_root(&mut metal_parent, other);
                if let Some(slot) = metal_parent.get_mut(other_root) {
                    *slot = first_root;
                }
            }
        }
    }
    let mut node_groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for &atom in &node_atoms {
        let parent = find_root(&mut metal_parent, atom);
        node_groups.entry(parent).or_default().push(atom);
    }
    let mut units: Vec<BuildingUnit> = Vec::new();
    for (_, mut atoms) in node_groups {
        atoms.sort_unstable();
        let elements = element_counts(&p1, &atoms)?;
        let node_metals: Vec<String> = {
            let mut set = BTreeSet::new();
            for &atom in &atoms {
                if let Some(site) = p1.sites().get(atom)
                    && element_info(site.element()).is_ok_and(|info| info.is_metal)
                {
                    set.insert(site.element().to_owned());
                }
            }
            set.into_iter().collect()
        };
        units.push(make_unit(&p1, UnitKind::Node, atoms, elements, node_metals));
    }
    // Linkers and free solvent.
    let mut assigned: BTreeSet<usize> = node_atoms.iter().copied().collect();
    for (index, group) in groups.iter().enumerate() {
        if index != main {
            let elements = element_counts(&p1, group)?;
            units.push(make_unit(
                &p1,
                UnitKind::Solvent,
                group.clone(),
                elements,
                Vec::new(),
            ));
            assigned.extend(group.iter().copied());
        }
    }
    let main_atoms: Vec<usize> = groups.get(main).cloned().unwrap_or_default();
    let remaining: Vec<usize> = main_atoms
        .iter()
        .copied()
        .filter(|index| !assigned.contains(index))
        .collect();
    let mut remaining_set: BTreeSet<usize> = remaining.iter().copied().collect();
    while let Some(&start) = remaining_set.iter().next() {
        let mut component = Vec::new();
        let mut stack = vec![start];
        remaining_set.remove(&start);
        while let Some(atom) = stack.pop() {
            component.push(atom);
            for &neighbour in adjacency.get(atom).map(Vec::as_slice).unwrap_or(&[]) {
                if remaining_set.remove(&neighbour) {
                    stack.push(neighbour);
                }
            }
        }
        component.sort_unstable();
        let elements = element_counts(&p1, &component)?;
        let is_linker = elements.contains_key("C")
            && component.iter().any(|atom| {
                adjacency
                    .get(*atom)
                    .map(Vec::as_slice)
                    .unwrap_or(&[])
                    .iter()
                    .any(|neighbour| node_atoms.contains(neighbour))
            });
        units.push(make_unit(
            &p1,
            if is_linker {
                UnitKind::Linker
            } else {
                UnitKind::Solvent
            },
            component,
            elements,
            Vec::new(),
        ));
    }
    // Connections.
    let mut atom_to_unit: BTreeMap<usize, usize> = BTreeMap::new();
    for (unit, entry) in units.iter().enumerate() {
        for &atom in &entry.atoms {
            atom_to_unit.insert(atom, unit);
        }
    }
    let mut connections: Vec<Vec<UnitConnection>> = vec![Vec::new(); units.len()];
    for bond in &graph.bonds {
        let (Some(&first_unit), Some(&second_unit)) =
            (atom_to_unit.get(&bond.i), atom_to_unit.get(&bond.j))
        else {
            continue;
        };
        if first_unit == second_unit || bond.i == bond.j {
            continue;
        }
        connections[first_unit].push(UnitConnection {
            atom: bond.i,
            other: bond.j,
            other_unit: second_unit,
            image: bond.image,
        });
        connections[second_unit].push(UnitConnection {
            atom: bond.j,
            other: bond.i,
            other_unit: first_unit,
            image: [-bond.image[0], -bond.image[1], -bond.image[2]],
        });
    }
    for (unit, list) in connections.into_iter().enumerate() {
        if let Some(entry) = units.get_mut(unit) {
            entry.connections = list;
            entry.span = connection_span(&p1, entry);
        }
    }
    let net = build_net(&p1, &units);
    let (topology, nearest_topologies) = identify_topology(&net);
    let framework_components = groups
        .iter()
        .filter(|group| {
            group.iter().any(|atom| {
                p1.sites().get(*atom).is_some_and(|site| {
                    element_info(site.element()).is_ok_and(|info| info.is_metal)
                })
            })
        })
        .count()
        .max(1);
    Ok(Deconstruction {
        structure: p1,
        units,
        net,
        topology,
        nearest_topologies,
        interpenetration: framework_components,
    })
}

fn make_unit(
    structure: &Structure,
    kind: UnitKind,
    atoms: Vec<usize>,
    elements: BTreeMap<String, usize>,
    metals: Vec<String>,
) -> BuildingUnit {
    let mut centroid = [0.0; 3];
    for &atom in &atoms {
        if let Some(site) = structure.sites().get(atom) {
            centroid[0] += site.frac()[0];
            centroid[1] += site.frac()[1];
            centroid[2] += site.frac()[2];
        }
    }
    let count = atoms.len().max(1) as f64;
    centroid = [
        centroid[0] / count,
        centroid[1] / count,
        centroid[2] / count,
    ];
    BuildingUnit {
        kind,
        formula: formula(&elements),
        atoms,
        elements,
        connections: Vec::new(),
        centroid,
        span: 0.0,
        metals,
    }
}

fn connection_span(structure: &Structure, unit: &BuildingUnit) -> f64 {
    let mut points = Vec::new();
    for connection in &unit.connections {
        if let Some(site) = structure.sites().get(connection.atom) {
            points.push(site.frac());
        }
    }
    let mut span = 0.0_f64;
    for (index, first) in points.iter().enumerate() {
        for second in points.iter().skip(index + 1) {
            span = span.max(structure.lattice().distance(*first, *second));
        }
    }
    span
}

fn build_net(structure: &Structure, units: &[BuildingUnit]) -> QuotientGraph {
    let mut vertices: Vec<usize> = Vec::new();
    for (index, unit) in units.iter().enumerate() {
        match unit.kind {
            UnitKind::Node => vertices.push(index),
            UnitKind::Linker if unit.connections.len() != 2 => vertices.push(index),
            _ => {}
        }
    }
    let vertex_of: BTreeMap<usize, usize> = vertices
        .iter()
        .enumerate()
        .map(|(vertex, unit)| (*unit, vertex))
        .collect();
    let mut nodes = Vec::with_capacity(vertices.len());
    let mut labels = Vec::with_capacity(vertices.len());
    for unit in &vertices {
        if let Some(entry) = units.get(*unit) {
            nodes.push(entry.centroid);
            labels.push(match entry.kind {
                UnitKind::Node => "M".to_owned(),
                _ => "L".to_owned(),
            });
        }
    }
    let mut edges = Vec::new();
    for unit in units {
        if unit.kind != UnitKind::Linker || unit.connections.len() != 2 {
            continue;
        }
        let mut targets = Vec::new();
        for connection in &unit.connections {
            if let Some(&vertex) = vertex_of.get(&connection.other_unit) {
                targets.push(vertex);
            }
        }
        targets.sort_unstable();
        targets.dedup();
        if targets.len() != 2 {
            continue;
        }
        let (Some(&a), Some(&b)) = (targets.first(), targets.get(1)) else {
            continue;
        };
        let (Some(first), Some(second)) = (nodes.get(a), nodes.get(b)) else {
            continue;
        };
        let delta = [
            second[0] - first[0],
            second[1] - first[1],
            second[2] - first[2],
        ];
        let shift = delta.map(|value| -value.round() as i32);
        let _ = structure;
        edges.push(NetEdge { a, b, shift });
    }
    edges.sort();
    edges.dedup();
    QuotientGraph {
        nodes,
        edges,
        labels,
    }
}

fn identify_topology(net: &QuotientGraph) -> (Vec<String>, Vec<String>) {
    if net.nodes.is_empty() || net.edges.is_empty() {
        return (Vec::new(), Vec::new());
    }
    match xtalforge_analysis::RcsrIndex::from_embedded() {
        Ok(index) => match index.identify(net) {
            Ok(identification) => (identification.names, identification.nearest),
            Err(_) => (Vec::new(), Vec::new()),
        },
        Err(_) => (Vec::new(), Vec::new()),
    }
}

fn find_root(parent: &mut [usize], mut node: usize) -> usize {
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
