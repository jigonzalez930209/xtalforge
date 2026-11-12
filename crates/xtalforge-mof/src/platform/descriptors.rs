//! Linker, node and framework descriptors (Task 10.2.5).

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use xtalforge_core::{Error, Structure, find_groups};

use crate::platform::deconstruct::{Deconstruction, UnitKind, deconstruct};

/// Descriptor of one linker.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LinkerDescriptor {
    /// Chemical formula.
    pub formula: String,
    /// Number of connections to nodes.
    pub connections: usize,
    /// Distance between the outer connection atoms in Angstrom.
    pub span: f64,
    /// Functional groups found on the linker.
    pub groups: BTreeMap<String, usize>,
}

/// Descriptor of one node.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NodeDescriptor {
    /// Chemical formula.
    pub formula: String,
    /// Metal elements.
    pub metals: Vec<String>,
    /// Number of metal atoms (nuclearity).
    pub nuclearity: usize,
    /// Number of connections to linkers.
    pub connectivity: usize,
}

/// The descriptor bundle used by MOF screening and ML.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FrameworkDescriptors {
    /// Framework formula.
    pub formula: String,
    /// Atoms in the P1 cell.
    pub n_atoms: usize,
    /// Cell volume in cubic Angstrom.
    pub cell_volume: f64,
    /// Density in g/cm3.
    pub density: f64,
    /// Distinct metal elements.
    pub metals: Vec<String>,
    /// RCSR topology names.
    pub topology: Vec<String>,
    /// Interpenetrating frameworks.
    pub interpenetration: usize,
    /// Node descriptors.
    pub nodes: Vec<NodeDescriptor>,
    /// Linker descriptors.
    pub linkers: Vec<LinkerDescriptor>,
    /// Functional-group counts over the whole framework.
    pub groups: BTreeMap<String, usize>,
}

/// Computes the linker, node and framework descriptors.
///
/// # Errors
///
/// Returns the errors of deconstruction, functional-group matching or the
/// density calculation.
pub fn framework_descriptors(structure: &Structure) -> Result<FrameworkDescriptors, Error> {
    let deconstruction = deconstruct(structure)?;
    descriptors_from(&deconstruction)
}

/// Computes descriptors from an existing deconstruction.
///
/// # Errors
///
/// Returns the errors of functional-group matching or density.
pub fn descriptors_from(deconstruction: &Deconstruction) -> Result<FrameworkDescriptors, Error> {
    let structure = &deconstruction.structure;
    let groups = find_groups(structure)?;
    let mut group_counts: BTreeMap<String, usize> = BTreeMap::new();
    for group in &groups {
        *group_counts
            .entry(group.kind.name().to_owned())
            .or_insert(0) += 1;
    }
    let mut nodes = Vec::new();
    let mut linkers = Vec::new();
    let mut metals: Vec<String> = Vec::new();
    for unit in &deconstruction.units {
        match unit.kind {
            UnitKind::Node => {
                for metal in &unit.metals {
                    if !metals.contains(metal) {
                        metals.push(metal.clone());
                    }
                }
                nodes.push(NodeDescriptor {
                    formula: unit.formula.clone(),
                    metals: unit.metals.clone(),
                    nuclearity: unit
                        .elements
                        .iter()
                        .filter(|(element, _)| unit.metals.contains(element))
                        .map(|(_, count)| *count)
                        .sum(),
                    connectivity: unit.connections.len(),
                });
            }
            UnitKind::Linker => {
                let mut unit_groups: BTreeMap<String, usize> = BTreeMap::new();
                for group in &groups {
                    if group.atoms.iter().all(|atom| unit.atoms.contains(atom)) {
                        *unit_groups.entry(group.kind.name().to_owned()).or_insert(0) += 1;
                    }
                }
                linkers.push(LinkerDescriptor {
                    formula: unit.formula.clone(),
                    connections: unit.connections.len(),
                    span: unit.span,
                    groups: unit_groups,
                });
            }
            UnitKind::Solvent => {}
        }
    }
    metals.sort();
    Ok(FrameworkDescriptors {
        formula: structure.formula(),
        n_atoms: structure.p1()?.n_atoms(),
        cell_volume: structure.lattice().volume(),
        density: structure.density()?,
        metals,
        topology: deconstruction.topology.clone(),
        interpenetration: deconstruction.interpenetration,
        nodes,
        linkers,
        groups: group_counts,
    })
}
