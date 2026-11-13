//! The "all-node" contraction rule of the MOF platform (Task 14.4.2).
//!
//! A metal and every non-carbon atom bonded to it, together with hydrogens on
//! those atoms, form one node. Metals that share such a bridge join the same
//! node. A carbon-containing unit bonded to exactly two nodes becomes one edge;
//! a unit bonded to three or more nodes is itself a node. Terminal groups and
//! free solvent are left out.
//!
//! See [`crate::contraction`] for the public entry point and the documented
//! centroid and shift conventions.

#![allow(clippy::indexing_slicing)] // bounded local node and edge arrays

use std::collections::{BTreeMap, BTreeSet};

use xtalforge_core::{BondGraph, Error, Structure, element_info};

use crate::contraction::{
    ContractedEdge, ContractedNet, ContractedNode, NodeKind, adjacency, centroid, find_root,
    finish, lattice_shift, make_edge, union,
};

pub(crate) fn contract_metal_clusters(
    structure: &Structure,
    bonds: &BondGraph,
) -> Result<ContractedNet, Error> {
    let n_atoms = structure.sites().len();
    let adjacency = adjacency(bonds, n_atoms);
    let mut is_carbon = vec![false; n_atoms];
    let mut is_hydrogen = vec![false; n_atoms];
    let mut metals = Vec::new();
    for (index, site) in structure.sites().iter().enumerate() {
        let info = element_info(site.element())?;
        is_carbon[index] = info.symbol == "C";
        is_hydrogen[index] = info.symbol == "H";
        if info.is_metal {
            metals.push(index);
        }
    }
    if metals.is_empty() {
        return Err(Error::InvalidInput(
            "the all-node rule needs at least one metal atom".to_owned(),
        ));
    }
    let mut parent: Vec<usize> = (0..n_atoms).collect();
    for &metal in &metals {
        for &neighbour in adjacency.get(metal).map(Vec::as_slice).unwrap_or(&[]) {
            if is_carbon.get(neighbour).copied().unwrap_or(false) {
                continue;
            }
            union(&mut parent, metal, neighbour);
            for &second in adjacency.get(neighbour).map(Vec::as_slice).unwrap_or(&[]) {
                if is_hydrogen.get(second).copied().unwrap_or(false) {
                    union(&mut parent, neighbour, second);
                }
            }
        }
    }
    let metal_roots: BTreeSet<usize> = metals
        .iter()
        .map(|metal| find_root(&mut parent, *metal))
        .collect();
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for atom in 0..n_atoms {
        let root = find_root(&mut parent, atom);
        if metal_roots.contains(&root) {
            groups.entry(root).or_default().push(atom);
        }
    }
    let mut nodes = Vec::new();
    let mut labels = Vec::new();
    let mut descriptions = Vec::new();
    let mut node_of_atom = vec![usize::MAX; n_atoms];
    for (_, mut atoms) in groups {
        atoms.sort_unstable();
        atoms.dedup();
        if atoms.is_empty() {
            continue;
        }
        let centroid = centroid(structure, &atoms)?;
        let vertex = nodes.len();
        for atom in &atoms {
            if let Some(slot) = node_of_atom.get_mut(*atom) {
                *slot = vertex;
            }
        }
        nodes.push(centroid);
        labels.push("M".to_owned());
        descriptions.push(ContractedNode {
            atoms,
            centroid,
            kind: NodeKind::MetalCluster,
            label: "M".to_owned(),
        });
    }
    if nodes.is_empty() {
        return Err(Error::InvalidInput(
            "the all-node rule produced no metal node".to_owned(),
        ));
    }
    let mut remaining: BTreeSet<usize> = (0..n_atoms)
        .filter(|index| node_of_atom[*index] == usize::MAX)
        .collect();
    let mut edges = Vec::new();
    while let Some(start) = remaining.iter().next().copied() {
        let mut component = Vec::new();
        let mut stack = vec![start];
        remaining.remove(&start);
        while let Some(atom) = stack.pop() {
            component.push(atom);
            for &neighbour in adjacency.get(atom).map(Vec::as_slice).unwrap_or(&[]) {
                if remaining.remove(&neighbour) {
                    stack.push(neighbour);
                }
            }
        }
        component.sort_unstable();
        let mut targets: BTreeSet<usize> = BTreeSet::new();
        for atom in &component {
            for &neighbour in adjacency.get(*atom).map(Vec::as_slice).unwrap_or(&[]) {
                if let Some(&vertex) = node_of_atom.get(neighbour)
                    && vertex != usize::MAX
                {
                    targets.insert(vertex);
                }
            }
        }
        if !component
            .iter()
            .any(|atom| is_carbon.get(*atom).copied().unwrap_or(false))
        {
            continue;
        }
        match targets.len() {
            2 => {
                let (Some(&a), Some(&b)) = (targets.first(), targets.last()) else {
                    continue;
                };
                if let (Some(first), Some(second)) = (nodes.get(a), nodes.get(b)) {
                    edges.push(ContractedEdge {
                        a,
                        b,
                        shift: lattice_shift(*first, *second),
                        atoms: component,
                    });
                }
            }
            3.. => {
                let centroid = centroid(structure, &component)?;
                let vertex = nodes.len();
                for atom in &component {
                    if let Some(slot) = node_of_atom.get_mut(*atom) {
                        *slot = vertex;
                    }
                }
                nodes.push(centroid);
                labels.push("L".to_owned());
                descriptions.push(ContractedNode {
                    atoms: component,
                    centroid,
                    kind: NodeKind::Linker,
                    label: "L".to_owned(),
                });
            }
            _ => {}
        }
    }
    let mut direct = Vec::new();
    for bond in &bonds.bonds {
        let (Some(&first), Some(&second)) = (node_of_atom.get(bond.i), node_of_atom.get(bond.j))
        else {
            continue;
        };
        if first == usize::MAX || second == usize::MAX || first == second {
            continue;
        }
        let (Some(a), Some(b)) = (nodes.get(first), nodes.get(second)) else {
            continue;
        };
        direct.push(make_edge(
            first,
            second,
            lattice_shift(*a, *b),
            vec![bond.i, bond.j],
        ));
    }
    direct.sort();
    direct.dedup();
    edges.extend(direct);
    finish(structure.clone(), nodes, labels, edges, descriptions, false)
}
