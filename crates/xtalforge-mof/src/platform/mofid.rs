//! MOFid and MOFkey identity strings (Task 10.2.3).

use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use xtalforge_build::sketch::{Sketch, SketchAtom};
use xtalforge_core::{Error, Structure};

use crate::platform::deconstruct::{Deconstruction, UnitKind, deconstruct};
use crate::platform::util::bond_adjacency;

/// The identity of a MOF.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MofId {
    /// RCSR topology name, or `unknown`.
    pub topology: String,
    /// Canonical SMILES of every distinct node.
    pub nodes: Vec<String>,
    /// Canonical SMILES of every distinct linker.
    pub linkers: Vec<String>,
    /// Canonical MOFid string.
    pub mofid: String,
    /// SHA-256 based MOFkey.
    pub mofkey: String,
    /// Number of interpenetrating frameworks.
    pub interpenetration: usize,
}

/// Builds the MOFid and MOFkey of a framework.
///
/// The string is a deterministic subset of the published MOFid syntax:
/// `MOFid:{topology}:{nodes}|{linkers}|{n}`. Nodes and linkers are sorted and
/// de-duplicated canonical graph strings (a depth-first SMILES writer over an
/// invariant-sorted graph).
///
/// # Errors
///
/// Returns the errors of deconstruction or of the SMILES writer.
pub fn mofid(structure: &Structure) -> Result<MofId, Error> {
    let deconstruction = deconstruct(structure)?;
    mofid_from_deconstruction(&deconstruction)
}

/// Builds the MOFid from an existing deconstruction.
///
/// # Errors
///
/// Returns the errors of the SMILES writer.
pub fn mofid_from_deconstruction(deconstruction: &Deconstruction) -> Result<MofId, Error> {
    let topology = deconstruction
        .topology
        .first()
        .cloned()
        .unwrap_or_else(|| "unknown".to_owned());
    let mut nodes = BTreeSet::new();
    let mut linkers = BTreeSet::new();
    for unit in &deconstruction.units {
        match unit.kind {
            UnitKind::Node => {
                nodes.insert(unit_smiles(&deconstruction.structure, unit)?);
            }
            UnitKind::Linker => {
                linkers.insert(unit_smiles(&deconstruction.structure, unit)?);
            }
            UnitKind::Solvent => {}
        }
    }
    let nodes: Vec<String> = nodes.into_iter().collect();
    let linkers: Vec<String> = linkers.into_iter().collect();
    let mofid = format!(
        "MOFid:{topology}:{}|{}|{}",
        nodes.join("."),
        linkers.join("."),
        deconstruction.interpenetration
    );
    let digest = Sha256::digest(mofid.as_bytes());
    let mofkey = digest
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect::<String>();
    Ok(MofId {
        topology,
        nodes,
        linkers,
        mofid,
        mofkey,
        interpenetration: deconstruction.interpenetration,
    })
}

/// Canonical SMILES of one unit, with external bonds cut.
///
/// # Errors
///
/// Returns the error of the SMILES writer.
pub fn unit_smiles(
    structure: &Structure,
    unit: &crate::platform::deconstruct::BuildingUnit,
) -> Result<String, Error> {
    let (p1, graph) = crate::platform::util::p1_and_bonds(structure)?;
    let adjacency = bond_adjacency(&graph, p1.sites().len());
    let index_of: BTreeMap<usize, usize> = unit
        .atoms
        .iter()
        .enumerate()
        .map(|(local, parent)| (*parent, local))
        .collect();
    let mut colours: Vec<(String, usize)> = unit
        .atoms
        .iter()
        .map(|parent| {
            let element = p1
                .sites()
                .get(*parent)
                .map(|site| site.element().to_owned())
                .unwrap_or_else(|| "X".to_owned());
            let degree = adjacency
                .get(*parent)
                .map(Vec::as_slice)
                .unwrap_or(&[])
                .iter()
                .filter(|neighbour| index_of.contains_key(neighbour))
                .count();
            (element, degree)
        })
        .collect();
    for _ in 0..5 {
        let snapshot = colours.clone();
        for (local, parent) in unit.atoms.iter().enumerate() {
            let mut neighbour_colours: Vec<(String, usize)> = adjacency
                .get(*parent)
                .map(Vec::as_slice)
                .unwrap_or(&[])
                .iter()
                .filter_map(|neighbour| index_of.get(neighbour).copied())
                .filter_map(|other| snapshot.get(other).cloned())
                .collect();
            neighbour_colours.sort();
            let mut colour = snapshot.get(local).cloned().unwrap_or_default();
            colour.0.push('|');
            for (element, degree) in neighbour_colours {
                colour.0.push_str(&format!("{element}:{degree},"));
            }
            if let Some(slot) = colours.get_mut(local) {
                *slot = colour;
            }
        }
    }
    let mut order: Vec<usize> = (0..unit.atoms.len()).collect();
    order.sort_by(|first, second| {
        colours
            .get(*first)
            .cmp(&colours.get(*second))
            .then_with(|| first.cmp(second))
    });
    let new_index: BTreeMap<usize, usize> = order
        .iter()
        .enumerate()
        .map(|(new, old)| (*old, new))
        .collect();
    let mut sketch = Sketch::new();
    for local in &order {
        let parent = unit.atoms.get(*local).copied().unwrap_or(0);
        let element = p1
            .sites()
            .get(parent)
            .map(|site| site.element().to_owned())
            .unwrap_or_else(|| "X".to_owned());
        sketch.add_atom(SketchAtom {
            element,
            position: [0.0, 0.0],
            charge: 0,
        });
    }
    let mut seen: BTreeSet<(usize, usize)> = BTreeSet::new();
    for (local, parent) in unit.atoms.iter().enumerate() {
        for &neighbour in adjacency.get(*parent).map(Vec::as_slice).unwrap_or(&[]) {
            let (Some(&here), Some(&there)) = (new_index.get(&local), index_of.get(&neighbour))
            else {
                continue;
            };
            let Some(&other) = new_index.get(&there) else {
                continue;
            };
            if here == other {
                continue;
            }
            let pair = if here < other {
                (here, other)
            } else {
                (other, here)
            };
            if !seen.insert(pair) {
                continue;
            }
            sketch.connect(here, other, 1.0)?;
        }
    }
    sketch.to_smiles()
}

/// Whether two MOFs have the same MOFkey.
///
/// # Errors
///
/// Returns the errors of [`mofid`].
pub fn same_mofkey(first: &Structure, second: &Structure) -> Result<bool, Error> {
    Ok(mofid(first)?.mofkey == mofid(second)?.mofkey)
}
