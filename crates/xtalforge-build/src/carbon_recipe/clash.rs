//! Removal of edge terminations that overlap another atom (Task 6.8.5).
//!
//! A terminating atom is placed along its carbon's free valence without
//! looking at the rest of the sheet, so across a narrow gap two terminations
//! can land on top of each other or on a carbon. Every termination group (an
//! atom bonded to the carbon graph plus the atoms bonded to it, such as the
//! hydrogen of a hydroxyl) is checked in placement order; a group with an atom
//! closer than [`CLASH`] to a non-bonded kept atom is dropped and its carbon
//! is left bare. Termination bonds are given their minimum lattice image.

#![allow(clippy::indexing_slicing)] // bounded site arrays

use xtalforge_core::{ExplicitBond, Lattice, Site};

/// Non-bonded atoms closer than this overlap, Angstrom.
pub(crate) const CLASH: f64 = 1.2;

/// Sites and bonds without the overlapping groups, and the dropped elements.
pub(crate) struct Pruned {
    pub(crate) sites: Vec<Site>,
    pub(crate) bonds: Vec<ExplicitBond>,
    /// Element symbol of every dropped atom.
    pub(crate) dropped: Vec<String>,
    /// Number of dropped groups, one per carbon left bare.
    pub(crate) groups: usize,
}

/// Drops the termination groups (sites from `first_terminal` on) that clash.
pub(crate) fn prune(
    lattice: &Lattice,
    sites: Vec<Site>,
    bonds: Vec<ExplicitBond>,
    first_terminal: usize,
) -> Pruned {
    let n = sites.len();
    // Root of every terminal atom: the terminal bonded to the carbon graph.
    let mut root: Vec<usize> = (0..n).collect();
    for bond in &bonds {
        if bond.j >= first_terminal && bond.j < n {
            root[bond.j] = if bond.i >= first_terminal {
                root[bond.i]
            } else {
                bond.j
            };
        }
    }
    let bonded = |a: usize, b: usize| {
        bonds
            .iter()
            .any(|bond| (bond.i == a && bond.j == b) || (bond.i == b && bond.j == a))
    };
    let mut keep = vec![true; n];
    for group in first_terminal..n {
        if root[group] != group {
            continue;
        }
        let members: Vec<usize> = (first_terminal..n)
            .filter(|atom| root[*atom] == group)
            .collect();
        let clashes = members.iter().any(|&atom| {
            (0..n).any(|other| {
                other != atom
                    && keep[other]
                    && root[other] != group
                    && !bonded(atom, other)
                    && (other < first_terminal || other < group)
                    && lattice.distance(sites[atom].frac(), sites[other].frac()) < CLASH
            })
        });
        if clashes {
            for atom in members {
                keep[atom] = false;
            }
        }
    }
    let mut remap = vec![usize::MAX; n];
    let mut kept_sites = Vec::with_capacity(n);
    let mut dropped = Vec::new();
    for (index, site) in sites.into_iter().enumerate() {
        if keep[index] {
            remap[index] = kept_sites.len();
            kept_sites.push(site);
        } else {
            dropped.push(site.element().to_owned());
        }
    }
    let groups = (first_terminal..n)
        .filter(|atom| root[*atom] == *atom && !keep[*atom])
        .count();
    let kept_bonds = bonds
        .into_iter()
        .filter_map(|bond| {
            let (i, j) = (*remap.get(bond.i)?, *remap.get(bond.j)?);
            if i == usize::MAX || j == usize::MAX {
                return None;
            }
            // A termination is placed in Cartesian space and may cross the
            // cell boundary; its bond takes the minimum image.
            let image = if bond.j >= first_terminal {
                let a = Lattice::wrap(kept_sites[i].frac());
                let b = Lattice::wrap(kept_sites[j].frac());
                [0, 1, 2].map(|axis| -((b[axis] - a[axis]).round()) as i32)
            } else {
                bond.image
            };
            Some(ExplicitBond {
                i,
                j,
                image,
                ..bond
            })
        })
        .collect();
    Pruned {
        sites: kept_sites,
        bonds: kept_bonds,
        dropped,
        groups,
    }
}
