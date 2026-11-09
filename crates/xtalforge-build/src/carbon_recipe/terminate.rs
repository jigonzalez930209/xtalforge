//! Edge termination and the carbon structure (Task 6.8.3).
//!
//! Edge carbons are the graph sites with fewer than three neighbours.
//! The recipe's ratios decide how many hydroxyls, carbonyls, fluorines
//! and hydrogens go on them; edge carbons beyond the requested count are
//! left bare. Bonds are stated, never perceived.

#![allow(clippy::indexing_slicing)] // bounded local arrays

use std::collections::BTreeMap;

use xtalforge_core::{Error, ExplicitBond, Lattice, Site, Structure, StructureBuilder};

use super::graph::CarbonGraph;
use super::net::{add, norm, sub};

/// Bond lengths used when a termination is placed, Angstrom.
const C_H: f64 = 1.09;
const C_F: f64 = 1.35;
const C_O: f64 = 1.36;
const O_H: f64 = 0.97;

/// The edge chemistry of a recipe, per carbon of the result.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Ratios {
    pub(crate) hydrogen: f64,
    pub(crate) fluorine: f64,
    pub(crate) oxygen: f64,
    pub(crate) ether: f64,
    pub(crate) hydroxyl: f64,
    pub(crate) carbonyl: f64,
}

impl Ratios {
    /// How many of each kind a sheet of `carbons` should carry.
    pub(crate) fn wanted(&self, carbons: usize) -> Wanted {
        let split = [self.ether, self.hydroxyl, self.carbonyl];
        let total: f64 = split.iter().sum();
        let fractions = if total > 0.0 {
            split.map(|value| value / total)
        } else {
            [0.0; 3]
        };
        let r_ether = self.oxygen * fractions[0];
        let ether = (r_ether * carbons as f64 / (1.0 + r_ether)).round() as usize;
        let carbons = carbons.saturating_sub(ether);
        let hydroxyl = (self.oxygen * fractions[1] * carbons as f64).round() as usize;
        let carbonyl = (self.oxygen * fractions[2] * carbons as f64).round() as usize;
        let fluorine = (self.fluorine * carbons as f64).round() as usize;
        let hydrogen = (self.hydrogen * carbons as f64).round() as usize;
        Wanted {
            ether,
            hydroxyl,
            carbonyl,
            fluorine,
            hydrogen: hydrogen.saturating_sub(hydroxyl),
        }
    }
}

/// Requested termination counts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Wanted {
    pub(crate) ether: usize,
    pub(crate) hydroxyl: usize,
    pub(crate) carbonyl: usize,
    pub(crate) fluorine: usize,
    pub(crate) hydrogen: usize,
}

/// A terminated graph.
pub(crate) struct Terminated {
    pub(crate) structure: Structure,
    pub(crate) carbon_structure: Structure,
    pub(crate) carbon_count: usize,
    pub(crate) hydrogen: usize,
    pub(crate) fluorine: usize,
    pub(crate) oxygen: usize,
    pub(crate) edge_carbons: usize,
    pub(crate) bare_edges: usize,
    pub(crate) pieces: usize,
    pub(crate) short: BTreeMap<String, usize>,
}

/// Terminates the edge carbons and builds the structure.
pub(crate) fn terminate(graph: &CarbonGraph, ratios: Ratios) -> Result<Terminated, Error> {
    let degree = graph.degree();
    let mut elements: Vec<String> = vec!["C".to_owned(); graph.len()];
    let mut edge_sites: Vec<usize> = (0..graph.len())
        .filter(|index| degree.get(*index).copied().unwrap_or(0) < 3)
        .collect();
    edge_sites.sort_by_key(|index| (degree.get(*index).copied().unwrap_or(0), *index));
    let edge_carbons = edge_sites.len();
    let wanted = ratios.wanted(graph.len());
    let mut cursor = 0usize;
    let mut used = vec![false; graph.len()];
    let mut take = |count: usize, used: &mut [bool]| -> Vec<usize> {
        let mut out = Vec::new();
        while out.len() < count {
            let Some(index) = edge_sites.get(cursor).copied() else {
                break;
            };
            cursor += 1;
            if used.get(index).copied().unwrap_or(true) {
                continue;
            }
            if let Some(slot) = used.get_mut(index) {
                *slot = true;
            }
            out.push(index);
        }
        out
    };
    let mut short = BTreeMap::new();
    let mut record = |kind: &str, wanted: usize, got: usize| {
        if got < wanted {
            short.insert(kind.to_owned(), wanted - got);
        }
    };
    let ether = take(wanted.ether, &mut used);
    record("ether", wanted.ether, ether.len());
    for index in &ether {
        if let Some(value) = elements.get_mut(*index) {
            *value = "O".to_owned();
        }
    }
    let hydroxyl = take(wanted.hydroxyl, &mut used);
    record("hydroxyl", wanted.hydroxyl, hydroxyl.len());
    let carbonyl = take(wanted.carbonyl, &mut used);
    record("carbonyl", wanted.carbonyl, carbonyl.len());
    let fluorine = take(wanted.fluorine, &mut used);
    record("fluorine", wanted.fluorine, fluorine.len());
    let hydrogen = take(wanted.hydrogen, &mut used);
    record("hydrogen", wanted.hydrogen, hydrogen.len());
    let bare_edges = edge_sites
        .len()
        .saturating_sub(used.iter().filter(|value| **value).count());
    let lattice = graph.lattice;
    let mut sites: Vec<Site> = Vec::with_capacity(graph.len() + 4 * edge_sites.len());
    for (index, center) in graph.centers.iter().enumerate() {
        let element = elements.get(index).map_or("C", String::as_str);
        sites.push(Site::new(format!("C{index}"), element, *center)?);
    }
    let mut bonds: Vec<ExplicitBond> = graph
        .bonds
        .iter()
        .map(|bond| ExplicitBond {
            i: bond.a,
            j: bond.b,
            image: bond.image,
            order: 1.0,
            remove: false,
        })
        .collect();
    let mut next = sites.len();
    let mut hydrogen_count = 0usize;
    let mut fluorine_count = 0usize;
    let mut oxygen_count = ether.len();
    for index in hydroxyl {
        let direction = outward(graph, &degree, index);
        let mut oxygen = sites[index].frac();
        let oxygen_cart = lattice.to_cart(oxygen);
        let position = add(oxygen_cart, direction.map(|value| value * C_O));
        oxygen = lattice.to_frac(position);
        sites.push(Site::new(format!("O{index}"), "O", oxygen)?);
        bonds.push(ExplicitBond {
            i: index,
            j: next,
            image: [0; 3],
            order: 1.0,
            remove: false,
        });
        let second = perpendicular(direction);
        let hydrogen_position = add(position, second.map(|value| value * O_H));
        sites.push(Site::new(
            format!("H{index}"),
            "H",
            lattice.to_frac(hydrogen_position),
        )?);
        bonds.push(ExplicitBond {
            i: next,
            j: next + 1,
            image: [0; 3],
            order: 1.0,
            remove: false,
        });
        next += 2;
        oxygen_count += 1;
        hydrogen_count += 1;
    }
    for index in carbonyl {
        let direction = outward(graph, &degree, index);
        let position = add(
            lattice.to_cart(sites[index].frac()),
            direction.map(|value| value * 1.23),
        );
        sites.push(Site::new(
            format!("O{index}"),
            "O",
            lattice.to_frac(position),
        )?);
        bonds.push(ExplicitBond {
            i: index,
            j: next,
            image: [0; 3],
            order: 2.0,
            remove: false,
        });
        next += 1;
        oxygen_count += 1;
    }
    for index in fluorine {
        let direction = outward(graph, &degree, index);
        let position = add(
            lattice.to_cart(sites[index].frac()),
            direction.map(|value| value * C_F),
        );
        sites.push(Site::new(
            format!("F{index}"),
            "F",
            lattice.to_frac(position),
        )?);
        bonds.push(ExplicitBond {
            i: index,
            j: next,
            image: [0; 3],
            order: 1.0,
            remove: false,
        });
        next += 1;
        fluorine_count += 1;
    }
    for index in hydrogen {
        let direction = outward(graph, &degree, index);
        let position = add(
            lattice.to_cart(sites[index].frac()),
            direction.map(|value| value * C_H),
        );
        sites.push(Site::new(
            format!("H{index}"),
            "H",
            lattice.to_frac(position),
        )?);
        bonds.push(ExplicitBond {
            i: index,
            j: next,
            image: [0; 3],
            order: 1.0,
            remove: false,
        });
        next += 1;
        hydrogen_count += 1;
    }
    let carbon_count = elements.iter().filter(|element| *element == "C").count();
    let structure = build(lattice, &sites, &bonds, "carbon")?;
    let carbon_sites: Vec<Site> = sites
        .iter()
        .take(graph.len())
        .filter(|site| site.element() == "C")
        .cloned()
        .collect();
    let carbon_index: Vec<usize> = sites
        .iter()
        .take(graph.len())
        .enumerate()
        .filter_map(|(index, site)| (site.element() == "C").then_some(index))
        .collect();
    let mut remap = vec![usize::MAX; graph.len()];
    for (new, old) in carbon_index.iter().enumerate() {
        if let Some(slot) = remap.get_mut(*old) {
            *slot = new;
        }
    }
    let carbon_bonds: Vec<ExplicitBond> = graph
        .bonds
        .iter()
        .filter_map(|bond| {
            let a = remap.get(bond.a).copied()?;
            let b = remap.get(bond.b).copied()?;
            if a == usize::MAX || b == usize::MAX {
                return None;
            }
            Some(ExplicitBond {
                i: a,
                j: b,
                image: bond.image,
                order: 1.0,
                remove: false,
            })
        })
        .collect();
    let carbon_structure = build(lattice, &carbon_sites, &carbon_bonds, "carbon-network")?;
    let pieces = graph
        .components()
        .iter()
        .collect::<std::collections::BTreeSet<_>>()
        .len();
    Ok(Terminated {
        structure,
        carbon_structure,
        carbon_count,
        hydrogen: hydrogen_count,
        fluorine: fluorine_count,
        oxygen: oxygen_count,
        edge_carbons,
        bare_edges,
        pieces,
        short,
    })
}

fn outward(graph: &CarbonGraph, degree: &[usize], index: usize) -> [f64; 3] {
    let matrix = graph.lattice.matrix();
    let mut sum = [0.0; 3];
    for bond in &graph.bonds {
        let other = if bond.a == index {
            Some((bond.b, bond.image))
        } else if bond.b == index {
            Some((bond.a, bond.image.map(|value| -value)))
        } else {
            None
        };
        if let Some((other, image)) = other {
            let start = graph.centers[index];
            let end = add(graph.centers[other], image.map(f64::from));
            let vector = transform_fraction(sub(end, start), matrix);
            sum = add(sum, vector);
        }
    }
    let normal = graph.normals.get(index).copied().unwrap_or([0.0; 3]);
    let direction = if degree.get(index).copied().unwrap_or(0) >= 2 && norm(normal) > 0.1 {
        normal
    } else {
        sum.map(|value| -value)
    };
    let length = norm(direction).max(1e-12);
    direction.map(|value| value / length)
}

fn perpendicular(direction: [f64; 3]) -> [f64; 3] {
    let reference = if direction[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let perpendicular = [
        direction[1] * reference[2] - direction[2] * reference[1],
        direction[2] * reference[0] - direction[0] * reference[2],
        direction[0] * reference[1] - direction[1] * reference[0],
    ];
    let length = norm(perpendicular).max(1e-12);
    perpendicular.map(|value| value / length)
}

fn transform_fraction(vector: [f64; 3], matrix: [[f64; 3]; 3]) -> [f64; 3] {
    [
        vector[0] * matrix[0][0] + vector[1] * matrix[1][0] + vector[2] * matrix[2][0],
        vector[0] * matrix[0][1] + vector[1] * matrix[1][1] + vector[2] * matrix[2][1],
        vector[0] * matrix[0][2] + vector[1] * matrix[1][2] + vector[2] * matrix[2][2],
    ]
}

fn build(
    lattice: Lattice,
    sites: &[Site],
    bonds: &[ExplicitBond],
    name: &str,
) -> Result<Structure, Error> {
    let mut builder = StructureBuilder::new(name, lattice);
    for site in sites {
        builder = builder.site(site.clone());
    }
    for bond in bonds {
        builder = builder.explicit_bond(*bond);
    }
    builder.build()
}
